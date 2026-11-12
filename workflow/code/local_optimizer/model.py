"""Ideal-sail RK4 segments and the multiple-shooting equations.

States use AU and AU/day; clock and cone use radians. The caller supplies the
number of holds and the time scale. Building the graph does not run Ipopt.
"""

from __future__ import annotations

import math

import casadi as ca
import numpy as np

AU_METRES = 149597870691.0
DAY_SECONDS = 86400.0
GM_SI = 1.32712440018e20
GM_AUD = GM_SI * DAY_SECONDS**2 / AU_METRES**3
SAIL_CUTOFF_AU = 5.0
DEPARTURE_MJD = 52654.5
LIGHTNESS = 0.09274704077328881


def build_segment(substeps):
    """Build the held-command flow, RHS and complete RK4-map Jacobian.

    The moving RTN frame is rebuilt at every RK stage. Lightness is a fixed
    model parameter; state, clock/cone and hold duration are the map inputs."""
    if isinstance(substeps, bool) or not isinstance(substeps, int):
        raise ValueError("substeps must be an integer, not a Boolean")
    if not 1 <= substeps <= 128:
        raise ValueError("substeps must be between 1 and 128 for this diagnostic")

    x = ca.SX.sym("x", 6)
    u = ca.SX.sym("u", 2)
    lightness = ca.SX.sym("lightness")
    r, v = x[:3], x[3:6]
    radius = ca.norm_2(r)
    radial = r / radius
    angular_momentum = ca.cross(r, v)
    normal_to_orbit = angular_momentum / ca.norm_2(angular_momentum)
    transverse = ca.cross(normal_to_orbit, radial)
    sail_normal = ca.cos(u[1]) * radial + ca.sin(u[1]) * (
        ca.cos(u[0]) * transverse + ca.sin(u[0]) * normal_to_orbit
    )
    factor = ca.if_else(radius > SAIL_CUTOFF_AU, 0.0, lightness * GM_AUD / radius**2)
    sail_acceleration = factor * ca.cos(u[1]) ** 2 * sail_normal
    derivative = ca.vertcat(v, -GM_AUD * r / radius**3 + sail_acceleration)
    rhs = ca.Function(
        f"ideal_sail_rhs_{substeps}",
        [x, u, lightness],
        [derivative],
        ["x", "u", "lightness"],
        ["dxdt"],
    )

    initial = ca.MX.sym("x", 6)
    control = ca.MX.sym("u", 2)
    duration = ca.MX.sym("h")
    light = ca.MX.sym("lightness")
    state = initial
    dt = duration / substeps
    for _ in range(substeps):
        # Rebuild the moving frame at each RK stage; freezing it changes the sail force.
        k1 = rhs(state, control, light)
        k2 = rhs(state + dt * k1 / 2, control, light)
        k3 = rhs(state + dt * k2 / 2, control, light)
        k4 = rhs(state + dt * k3, control, light)
        state = state + dt * (k1 + 2 * k2 + 2 * k3 + k4) / 6

    inputs = [initial, control, duration, light]
    names = ["x", "u", "h", "lightness"]
    flow = ca.Function(f"ideal_sail_hold_{substeps}", inputs, [state], names, ["x_end"])
    jacobian = ca.Function(
        f"ideal_sail_hold_jacobian_{substeps}",
        inputs,
        [
            ca.horzcat(
                ca.jacobian(state, initial),
                ca.jacobian(state, control),
                ca.jacobian(state, duration),
            )
        ],
        names,
        ["d_x_end_d_x_u_h"],
    )
    return {"flow": flow, "jacobian": jacobian, "rhs": rhs}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def _scalar(value, name):
    require(
        isinstance(value, (int, float, np.integer, np.floating))
        and not isinstance(value, (bool, np.bool_)),
        name + " must be numeric",
    )
    result = float(value)
    require(math.isfinite(result), name + " must be finite")
    return result


def _array(value, shape, name):
    result = np.asarray(value, dtype=float)
    require(result.shape == shape and np.isfinite(result).all(), "Wrong shape/nonfinite " + name)
    return result.copy()


def _segment():
    return build_segment(16)


def checked_discretisation(settings):
    """Validate the mesh declaration without resizing or retiming seed data."""
    require(isinstance(settings, dict), "Settings must be a contract dictionary")
    n = settings.get("holds")
    require(
        isinstance(n, (int, np.integer)) and not isinstance(n, (bool, np.bool_)) and n > 0,
        "holds must be a positive integer",
    )
    n = int(n)
    state_count = settings.get("state_count", n + 1)
    require(
        isinstance(state_count, (int, np.integer))
        and not isinstance(state_count, (bool, np.bool_))
        and state_count == n + 1,
        "state_count must equal holds + 1",
    )
    duration_scale = _scalar(settings.get("duration_scale_days"), "duration_scale_days")
    require(duration_scale > 0.0, "duration_scale_days must be positive")
    return n, duration_scale


def build_problem(
    initial_state, seed_states, seed_controls, lightness, arrival_model, settings, seed_duration
):
    """Build the shooting graph, bounds and starting vector without solving.

    The supplied history is not retimed here. The departure is fixed, while
    later nodes and held angles are variables. Node guards do not replace
    the separate continuous-flight checks."""
    require(ca.__version__ == "3.7.2", "CasADi version changed")
    n, duration_scale = checked_discretisation(settings)
    fixed_numeric = {
        "rk_substeps": 16,
        "departure_MJD": DEPARTURE_MJD,
        "position_limit_km": 100000.0,
        "velocity_limit_km_s": 0.1,
        "position_reserve_km": 1.0,
        "velocity_reserve_m_s": 0.01,
        "velocity_box_AU_per_day": 0.1,
        "minimum_angular_momentum_AU2_day": 1e-8,
        "angular_momentum_scale_AU2_day": 0.02,
        "minimum_cos_latitude": 1e-8,
        "iterations": 300,
        "solver_wall_seconds": 90.0,
        "solver_cpu_seconds": 90.0,
        "tol": 1e-8,
        "constr_viol_tol": 1e-9,
        "acceptable_iter": 0,
        "bound_relax_factor": 0.0,
        "mumps_pivot_order": 0,
    }
    for key, expected in fixed_numeric.items():
        require(
            key in settings and _scalar(settings[key], key) == expected,
            "Frozen continuation setting changed: " + key,
        )
    fixed_strings = {
        "hessian_approximation": "exact",
        "linear_solver": "mumps",
        "honor_original_bounds": "yes",
        "nlp_scaling_method": "gradient-based",
    }
    for key, expected in fixed_strings.items():
        require(settings.get(key) == expected, "Frozen continuation setting changed: " + key)
    scale = _array(settings.get("state_scale"), (6,), "state scale")
    require(np.array_equal(scale, [1.0, 1.0, 1.0, 0.02, 0.02, 0.02]), "State scaling changed")
    radii = _array(settings.get("radius_bounds_AU"), (2,), "radius bounds")
    require(np.array_equal(radii, [0.2, 1.5]), "Original radial guards changed")
    duration_bounds = _array(settings.get("duration_bounds_days"), (2,), "duration bounds")
    require(
        0.5 <= duration_bounds[0] < duration_bounds[1] <= 600.0,
        "Duration lies outside the declared development domain",
    )
    require(
        isinstance(arrival_model, dict)
        and np.array_equal(duration_bounds, arrival_model.get("duration_bounds_days")),
        "Solver bounds must equal the selected target model bounds",
    )
    cell = _array(settings.get("arrival_cell_MJD"), (2,), "arrival cell")
    provenance = arrival_model.get("provenance", {})
    require(
        np.array_equal(cell, arrival_model.get("cell"))
        and np.array_equal(cell, provenance.get("cell_MJD")),
        "Arrival coefficient cell mismatch",
    )
    require(
        provenance.get("departure_MJD") == DEPARTURE_MJD
        and provenance.get("target") == "mercury"
        and provenance.get("global_dispatch_differentiated") is False
        and provenance.get("symbolic_bounds_enforced_by_caller") is True,
        "Unmatched exact-cell target model",
    )
    require(
        cell[0]
        <= DEPARTURE_MJD + duration_bounds[0]
        < DEPARTURE_MJD + duration_bounds[1]
        <= cell[1],
        "NLP duration window crosses target coefficient cells",
    )
    lightness = _scalar(lightness, "Lightness")
    require(lightness == LIGHTNESS, "Native-exported lightness changed; no free throttle")
    seed_duration = _scalar(seed_duration, "Seed duration")
    require(
        duration_bounds[0] <= seed_duration <= duration_bounds[1],
        "Seed duration is outside the selected target cell; no silent clipping",
    )
    initial = _array(initial_state, (6,), "initial state")
    states = _array(seed_states, (n + 1, 6), "seed states")
    controls = _array(seed_controls, (n, 2), "seed controls")
    require(np.array_equal(states[0], initial), "Original initial state must remain fixed")
    require(
        np.all(np.abs(controls[:, 0]) <= math.pi)
        and np.all((controls[:, 1] >= 0.0) & (controls[:, 1] <= math.pi / 2)),
        "Seed controls violate physical angle bounds",
    )
    require(
        np.all(np.abs(states[:, :3]) <= radii[1])
        and np.all(np.abs(states[:, 3:]) <= settings["velocity_box_AU_per_day"]),
        "Seed violates the Cartesian diagnostic box",
    )

    segment = _segment()
    w = ca.MX.sym("time_variables", 8 * n + 1)
    scaled_nodes = ca.reshape(w[: 6 * n], 6, n)
    commands = ca.reshape(w[6 * n : 8 * n], 2, n)
    total = w[-1] * duration_scale
    nodes = ca.times(ca.repmat(ca.DM(scale), 1, n), scaled_nodes)
    continuity = []
    for i in range(n):
        start = ca.DM(initial) if i == 0 else nodes[:, i - 1]
        # Nodes may move during the solve, but this equality must join the pieces.
        predicted = segment["flow"](start, commands[:, i], total / n, lightness)
        # Each piece must meet the next. No borrowed kilometres.
        continuity.append(scaled_nodes[:, i] - predicted / ca.DM(scale))
    continuity = ca.vertcat(*continuity)
    # Changing flight time moves Mercury too; do not keep the old arrival target.
    target = arrival_model["state"](total)
    endpoint = nodes[:, -1] - target
    position_nlp_km = settings["position_limit_km"] - settings["position_reserve_km"]
    velocity_nlp_km_s = settings["velocity_limit_km_s"] - settings["velocity_reserve_m_s"] / 1000.0
    require(position_nlp_km > 0.0 and velocity_nlp_km_s > 0.0, "Invalid reserve")
    optimisation_scale = ca.DM(
        [position_nlp_km * 1000 / AU_METRES] * 3
        + [velocity_nlp_km_s * 1000 * DAY_SECONDS / AU_METRES] * 3
    )
    acceptance_scale = ca.DM(
        [settings["position_limit_km"] * 1000 / AU_METRES] * 3
        + [settings["velocity_limit_km_s"] * 1000 * DAY_SECONDS / AU_METRES] * 3
    )
    endpoint_scaled = endpoint / optimisation_scale
    encounter = ca.vertcat(
        ca.dot(endpoint_scaled[:3], endpoint_scaled[:3]),
        ca.dot(endpoint_scaled[3:], endpoint_scaled[3:]),
    )
    objective = total / duration_scale

    # These constraints apply to nodes only. Independent RK-stage/path checks
    # remain mandatory; no clipping/regularisation hides singular trial states.
    radius_squared, angular_momentum_ratio, cosine_latitude = [], [], []
    for i in range(n):
        r, v = nodes[:3, i], nodes[3:6, i]
        radius_squared.append(ca.dot(r, r))
        angular_momentum_ratio.append(
            ca.norm_2(ca.cross(r, v)) / settings["angular_momentum_scale_AU2_day"]
        )
        cosine_latitude.append(ca.norm_2(r[:2]) / ca.norm_2(r))
    domain = ca.vertcat(*radius_squared, *angular_momentum_ratio, *cosine_latitude)
    constraints = ca.vertcat(continuity, domain, encounter)
    w0 = np.r_[(states[1:] / scale).ravel(), controls.ravel(), seed_duration / duration_scale]
    node_box = np.r_[[radii[1]] * 3, [settings["velocity_box_AU_per_day"]] * 3] / scale
    lbx = np.r_[
        np.tile(-node_box, n), np.tile([-math.pi, 0.0], n), duration_bounds[0] / duration_scale
    ]
    ubx = np.r_[
        np.tile(node_box, n),
        np.tile([math.pi, math.pi / 2], n),
        duration_bounds[1] / duration_scale,
    ]
    lbg = np.r_[
        np.zeros(6 * n),
        np.full(n, radii[0] ** 2),
        np.full(
            n,
            settings["minimum_angular_momentum_AU2_day"]
            / settings["angular_momentum_scale_AU2_day"],
        ),
        np.full(n, settings["minimum_cos_latitude"]),
        [-np.inf, -np.inf],
    ]
    # The former cosine upper bound was identically <= 1 on its nonsingular
    # domain. Remove that redundant inequality, not its lower frame safeguard.
    ubg = np.r_[np.zeros(6 * n), np.full(n, radii[1] ** 2), np.full(2 * n, np.inf), [1.0, 1.0]]
    require(np.all(w0 >= lbx) and np.all(w0 <= ubx), "Seed violates original variable bounds")
    functions = {
        "f": ca.Function("time_f", [w], [objective]),
        "g": ca.Function("time_g", [w], [constraints]),
        "gradient_f": ca.Function("time_gradient_f", [w], [ca.gradient(objective, w)]),
        "jac_g": ca.Function("time_jac_g", [w], [ca.jacobian(constraints, w)]),
        "continuity": ca.Function("time_continuity", [w], [continuity]),
        "endpoint_scaled": ca.Function("time_endpoint_scaled", [w], [endpoint_scaled]),
        "endpoint_acceptance_scaled": ca.Function(
            "time_endpoint_acceptance_scaled", [w], [endpoint / acceptance_scale]
        ),
        "endpoint_norms_squared": ca.Function("time_encounter", [w], [encounter]),
        "endpoint_physical": ca.Function("time_endpoint_physical", [w], [endpoint]),
        "nodes_physical": ca.Function("time_nodes_physical", [w], [nodes]),
    }
    solver_options = {
        "print_time": True,
        "error_on_fail": False,
        "ipopt.print_level": 5,
        "ipopt.max_iter": 300,
        "ipopt.max_wall_time": 90.0,
        "ipopt.max_cpu_time": 90.0,
        "ipopt.tol": 1e-8,
        "ipopt.constr_viol_tol": 1e-9,
        "ipopt.acceptable_iter": 0,
        "ipopt.hessian_approximation": "exact",
        "ipopt.linear_solver": "mumps",
        "ipopt.mumps_pivot_order": 0,
        "ipopt.bound_relax_factor": 0.0,
        "ipopt.honor_original_bounds": "yes",
        "ipopt.nlp_scaling_method": "gradient-based",
        "ipopt.check_derivatives_for_naninf": "yes",
    }
    return {
        "nlp": {"x": w, "f": objective, "g": constraints},
        "w0": w0,
        "lbx": lbx,
        "ubx": ubx,
        "lbg": lbg,
        "ubg": ubg,
        "functions": functions,
        "solver_options": solver_options,
        "layout": {
            "holds": n,
            "state_count": n + 1,
            "variables": 8 * n + 1,
            "constraints": 9 * n + 2,
            "state_scale": scale.tolist(),
            "duration_scale_days": duration_scale,
            "nodes": [0, 6 * n],
            "controls": [6 * n, 8 * n],
            "duration": 8 * n,
            "continuity": [0, 6 * n],
            "radius_squared": [6 * n, 7 * n],
            "angular_momentum_norm_scaled": [7 * n, 8 * n],
            "cos_latitude": [8 * n, 9 * n],
            "encounter_norms_squared": [9 * n, 9 * n + 2],
        },
        "model_provenance": {
            "segment_source": "code/local_optimizer/model.py",
            "arrival": dict(provenance),
        },
        "numerical_reserve": {
            "position_km": settings["position_reserve_km"],
            "velocity_m_s": settings["velocity_reserve_m_s"],
            "NLP_position_limit_km": position_nlp_km,
            "NLP_velocity_limit_km_s": velocity_nlp_km_s,
            "mission_position_limit_km": settings["position_limit_km"],
            "mission_velocity_limit_km_s": settings["velocity_limit_km_s"],
            "mission_acceptance_limits_unchanged": True,
        },
        "limitations": {
            "node_guards_only": True,
            "no_trial_domain_guarantee": True,
            "window_restricted_local_time_minimisation": True,
            "global_minimum_claim": False,
            "encounter_is_hard_constraint": True,
            "redundant_cos_latitude_upper_bound_removed": True,
            "continuity_nodes_not_a_propagated_flight": True,
            "native_replay_required_for_mission_acceptance": True,
        },
    }


def unpack_solution(vector, state_scale, *, holds, duration_scale_days):
    """Convert the declared scaled vector back to physical nodes, angles and time.

    Callers pass the hold count and time scale from the problem layout.
    This conversion does not establish feasibility."""
    n, duration_scale = checked_discretisation(
        {"holds": holds, "duration_scale_days": duration_scale_days}
    )
    point = _array(vector, (8 * n + 1,), "solution vector")
    scale = _array(state_scale, (6,), "state scale")
    require(np.array_equal(scale, [1.0, 1.0, 1.0, 0.02, 0.02, 0.02]), "State scaling changed")
    return (
        point[: 6 * n].reshape(n, 6) * scale,
        point[6 * n : 8 * n].reshape(n, 2),
        float(point[-1] * duration_scale),
    )
