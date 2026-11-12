"""Build and solve the 600-hold Mercury restoration and time problems.

The workflow owns native calls and run files. This module supplies the local
NLP, derivatives, one Ipopt call and the independent candidate checks.
"""

from __future__ import annotations

import copy
import hashlib
import json
import math
from pathlib import Path
import time

import casadi as ca
import numpy as np

from . import controls, ephemeris, model, objective, records, refinement, transfer, validation

SCHEMA = "fresh600_local_v5"
HOLDS = 600
HORIZON = 600.0
CELL = (53248.0, 53256.0)
DURATION = (593.5, 600.0)
STATE_SCALE = np.array([1.0, 1.0, 1.0, 0.02, 0.02, 0.02])
PHASES = ("restore", "quality")


def require(ok, message):
    if not ok:
        raise ValueError(message)


def _profile(profile=None):
    """This adapter freezes Mercury physics, although the objective is portable."""
    p = objective.validate_profile(profile)
    original = objective.default_profile()
    for key in (
        "reference_radius_m",
        "solar_mu_m3_s2",
        "global_position_limit_m",
        "global_velocity_limit_m_s",
        "numerical_position_gate_m",
        "numerical_velocity_gate_m_s",
    ):
        require(p[key] == original[key], "Mercury physics/global/numerical setting changed: " + key)
    return p


def protocol(mission_profile=None, continuation_domain=None, control_policy=None):
    """Versioned local development goals are separate from global ENC limits."""
    p = _profile(mission_profile)
    domains = ephemeris.campaign_domains(continuation_domain)
    domain_contract = copy.deepcopy(domains[0]["continuation_domain"])
    checked_controls = controls.validate_policy(control_policy)
    require(
        domain_contract is None or domain_contract["maximum_duration_days"] == HORIZON,
        "Fresh600 continuation must retain the full600-day export horizon",
    )
    long_campaign = domain_contract is not None
    attempts = 3 if long_campaign else 1
    budget = 900.0 if long_campaign else 600.0
    return {
        "schema": SCHEMA,
        "holds": HOLDS,
        "state_count": HOLDS + 1,
        "horizon_days": HORIZON,
        "arrival_cell_MJD": list(CELL),
        "duration_bounds_days": list(DURATION),
        "phases": list(PHASES),
        "continuation_domain": domain_contract,
        "control_policy": checked_controls,
        "arrival_domains": domains,
        "maximum_quality_solver_calls": len(domains) * attempts,
        "solver_calls_per_cell": attempts,
        "cell_capacity_is_not_promised_execution_count": True,
        "boundary_handover_max_retime_days": 1e-6,
        "mission_profile": p,
        "restore_objective": "0.5*sum((endpoint_error/restoration_conditioning_scale)^2)",
        "quality_objective": (
            "T_days/t_dyn_days; t_dyn_days=sqrt(r0^3/mu_sun)/86400; "
            "hard local-quality constraints"
        ),
        "time_scale_days": objective.dynamical_time_days(p),
        "native_fitness": "unchanged; reporting metric, not this local objective",
        "solver_calls_per_invocation": 1,
        "maximum_solver_calls": 1 + len(domains) * attempts,
        "iterations_per_phase": 3000 if long_campaign else 2000,
        "solver_wall_seconds_per_phase": budget,
        "solver_cpu_seconds_per_phase": budget,
        "per_call_solver_budget_seconds_range": [60.0, budget],
        "tol": 1e-10,
        "constr_viol_tol": 1e-12,
        "acceptable_iter": 0,
        "admission_policy": refinement.ADMISSION_POLICY,
        "budget_admission_scope": (
            "quality phase only; original bounds, objective, constraints "
            "and continuous validation required"
        ),
        "restoration_host_acceptance_required_before_quality": True,
        "restoration_native_replay": False,
        "final_native_replay_required": True,
        # Only continue a local result that passes the checks; do not restart ENC.
        "automatic_retry": long_campaign,
        "cell_hopping": True,
        "retry_scope": (
            "declared budget-feasible progress only, bounded by campaign clock"
            if long_campaign
            else "none"
        ),
        "cell_handover_is_explicit_and_host_validated": True,
        "fixed_restoration_duration_days": HORIZON,
        "position_limit_km": p["position_target_m"] / 1000.0,
        "velocity_limit_m_s": p["velocity_target_m_s"],
        "position_reserve_km": p["position_nlp_reserve_m"] / 1000.0,
        "velocity_reserve_m_s": p["velocity_nlp_reserve_m_s"],
        "restoration_scales_are_not_final_quality_limits": True,
        "whole_path_domain": "sampled physical and numerical gates, not an integrated running cost",
        "global_optimum_claim": False,
    }


def _runtime():
    require(
        ca.__version__ == "3.7.2" and np.__version__ == "2.5.3",
        "Reviewed numerical runtime changed",
    )


def _digest(value):
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()
    ).hexdigest()


def _identity(candidate, mission_profile=None):
    p = _profile(candidate.get("mission_profile") if mission_profile is None else mission_profile)
    require(candidate.get("mission_profile") == p, "Candidate mission profile differs")
    require(
        type(candidate.get("run_id")) is str
        and candidate["run_id"]
        and all(c.isalnum() or c in "-_" for c in candidate["run_id"]),
        "Invalid run identity",
    )
    digest = candidate.get("origin_controller_sha256")
    require(
        type(digest) is str and len(digest) == 64 and all(c in "0123456789abcdef" for c in digest),
        "Missing controller digest",
    )
    domain_contract = ephemeris.validate_domain(candidate.get("continuation_domain"))
    control_policy = controls.validate_policy(candidate.get("control_policy"))
    require(
        candidate.get("local_protocol_sha256")
        == _digest(protocol(p, domain_contract, control_policy)),
        "Local protocol changed",
    )
    return p


def load_full_export(
    export_dir,
    ephemeris_path,
    *,
    controller_sha256,
    run_id,
    mission_profile=None,
    continuation_domain=None,
    control_policy=None,
):
    """Check all600holds; no padding, truncation or ordinary-CSV tail synthesis.

    The controller digest is a required lineage label, not proof of fresh search.
    The parent must bind it to its same-run native search/export receipts.
    """
    profile = _profile(mission_profile)
    domain_contract = ephemeris.validate_domain(continuation_domain)
    control_policy = controls.validate_policy(control_policy)
    configuration = protocol(profile, domain_contract, control_policy)
    directory = Path(export_dir).resolve(strict=True)
    ephemeris_file = Path(ephemeris_path).resolve(strict=True)
    names = ("bridge_initial.csv", "bridge_states.csv", "bridge_commands.csv")
    paths = [directory / name for name in names] + [ephemeris_file]
    require(
        all(path.is_file() and not path.is_symlink() for path in paths),
        "Missing/linked full-export input",
    )
    bindings = {str(path): records.sha(path) for path in paths}
    initial = transfer._read_numeric(paths[0], transfer.INITIAL_FIELDS, 1)
    require(
        len(initial) == 1 and list(initial[0].values()) == transfer.FIXED_INITIAL,
        "Original initial state/model changed",
    )
    states = transfer._read_numeric(paths[1], transfer.STATE_FIELDS, HOLDS + 1)
    commands = transfer._read_numeric(paths[2], transfer.COMMAND_FIELDS, HOLDS + 1)
    require(len(states) == len(commands) == HOLDS + 1, "Full600export requires601boundaries")
    exposure = 0.0
    for i, (state, command) in enumerate(zip(states, commands)):
        require(
            state["command_index"] == i
            and state["MJD_d"] == command["MJD_d"] == transfer.DEPARTURE_MJD + i,
            "Full-export clock/index differs",
        )
        require(
            all(state[key] == command[key] for key in transfer.COMMAND_FIELDS),
            "Full-export control/state mismatch",
        )
        require(
            state["lightness"] == transfer.LIGHTNESS
            and state["ideal_sail_flag"] == 1.0
            and state["mp_kg"] == 0.0
            and state["srd_native"] >= exposure,
            "Full-export model or auxiliary state differs",
        )
        transfer._angles(command["clock_rad"], command["cone_rad"])
        transfer._check_polar(state)
        exposure = state["srd_native"]
    require(
        [states[0][key] for key in transfer.INITIAL_FIELDS] == transfer.FIXED_INITIAL,
        "First native boundary differs",
    )
    # Check target coverage at every actual epoch without fabricating target rows.
    eph = ephemeris.Ephemeris(ephemeris_file)
    for command in commands:
        require(
            np.isfinite(eph.heliocentric("mercury", command["MJD_d"])).all(), "Target coverage gap"
        )
    seed = {
        "schema": SCHEMA,
        "run_id": run_id,
        "origin_controller_sha256": controller_sha256,
        "local_protocol_sha256": _digest(configuration),
        "mission_profile": profile,
        "continuation_domain": domain_contract,
        "control_policy": control_policy,
        "departure_MJD": transfer.DEPARTURE_MJD,
        "duration_days": HORIZON,
        "lightness": transfer.LIGHTNESS,
        "holds": HOLDS,
        "state_count": HOLDS + 1,
        "nodes_cartesian": [[row[key] for key in transfer.CARTESIAN_FIELDS] for row in states],
        "controls_clock_cone_rad": [[row["clock_rad"], row["cone_rad"]] for row in commands[:-1]],
        "terminal_command_clock_cone_rad": [commands[-1]["clock_rad"], commands[-1]["cone_rad"]],
        "ephemeris_path": str(ephemeris_file),
        "source_sha256": bindings,
        "source_export_directory": str(directory),
        "source_duration_days": HORIZON,
        "discarded_tail_intervals": 0,
        "same_run_lineage_checked_by_adapter": False,
        "physical_transfer_only": True,
        "native_mission_accepted": False,
    }
    _identity(seed)
    require(
        all(records.sha(path) == digest for path, digest in bindings.items()),
        "Export changed while reading",
    )
    return seed


# Mercury does not wait at the old arrival date while we change T.
def build_arrival(ephemeris_path, cell_MJD=None, continuation_domain=None):
    """One explicitly selected allowed cell; no symbolic cross-cell dispatch."""
    _runtime()
    descriptor = ephemeris.descriptor(ephemeris_path, cell_MJD, continuation_domain)
    return ephemeris.build_cell(ephemeris_path, descriptor)


def baseline_candidate(seed, ephemeris_path, *, mission_profile=None):
    """Return the unchanged full history, not an accepted trajectory."""
    profile = _identity(seed, mission_profile)
    require(
        seed.get("schema") == SCHEMA
        and seed.get("holds") == HOLDS
        and seed.get("state_count") == HOLDS + 1
        and seed.get("duration_days") == HORIZON,
        "Wrong full-export seed",
    )
    require(
        seed.get("ephemeris_path") == str(Path(ephemeris_path).resolve()),
        "Seed ephemeris path differs",
    )
    bindings = seed.get("source_sha256", {})
    require(
        bool(bindings) and all(records.sha(path) == digest for path, digest in bindings.items()),
        "Source export/ephemeris changed",
    )
    rebuilt = load_full_export(
        seed["source_export_directory"],
        ephemeris_path,
        controller_sha256=seed["origin_controller_sha256"],
        run_id=seed["run_id"],
        mission_profile=profile,
        continuation_domain=seed.get("continuation_domain"),
        control_policy=seed.get("control_policy"),
    )
    require(seed == rebuilt, "Seed differs from actual full export; no silent array replacement")
    arrival = build_arrival(ephemeris_path, continuation_domain=seed.get("continuation_domain"))
    candidate = {
        key: copy.deepcopy(seed[key])
        for key in (
            "run_id",
            "origin_controller_sha256",
            "local_protocol_sha256",
            "mission_profile",
            "continuation_domain",
            "control_policy",
            "departure_MJD",
            "duration_days",
            "lightness",
            "holds",
            "state_count",
            "nodes_cartesian",
            "controls_clock_cone_rad",
            "terminal_command_clock_cone_rad",
        )
    }
    candidate.update(
        schema=1,
        arrival_cell_MJD=list(CELL),
        duration_bounds_days=list(DURATION),
        solver_target_cartesian=np.asarray(arrival["evaluate"](HORIZON)).reshape(6).tolist(),
        solver_phase="baseline",
        solver_status="NOT_RUN",
        native_mission_accepted=False,
    )
    validation.checked_candidate(candidate, "baseline")
    return candidate


def with_precise_nodes(candidate, info):
    """Use the already audited precise baseline nodes; commands remain unchanged."""
    _identity(candidate)
    validation.checked_candidate(candidate, "baseline")
    require(
        type(info) is dict
        and info.get("scientific_records_bound") is True
        and info.get("stage") == "baseline"
        and info.get("candidate") == candidate
        and info.get("commands") == validation.expected_commands(candidate),
        "Precise baseline is not bound to this candidate",
    )
    nodes = validation._array(info.get("native_cartesian_states"), (HOLDS + 1, 6), "precise nodes")
    require(np.array_equal(nodes[0], validation.INITIAL[:6]), "Precise initial state differs")
    result = copy.deepcopy(candidate)
    result["nodes_cartesian"] = nodes.tolist()
    result["starting_nodes"] = "same-run precise native baseline"
    return result


def _settings(cell_MJD=None, continuation_domain=None):
    # Start with the shooting equations, then set the phase objective, arrival
    # caps and solver budget explicitly in build_problem().
    selected = ephemeris.domain(cell_MJD, continuation_domain)
    return {
        "holds": HOLDS,
        "state_count": HOLDS + 1,
        "duration_scale_days": HORIZON,
        "arrival_cell_MJD": selected["cell_MJD"],
        "duration_bounds_days": selected["duration_bounds_days"],
        "rk_substeps": 16,
        "departure_MJD": transfer.DEPARTURE_MJD,
        "position_limit_km": 100000.0,
        "velocity_limit_km_s": 0.1,
        "position_reserve_km": 1.0,
        "velocity_reserve_m_s": 0.01,
        "state_scale": STATE_SCALE.tolist(),
        "radius_bounds_AU": [0.2, 1.5],
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
        "hessian_approximation": "exact",
        "linear_solver": "mumps",
        "mumps_pivot_order": 0,
        "bound_relax_factor": 0.0,
        "honor_original_bounds": "yes",
        "nlp_scaling_method": "gradient-based",
    }


def _solver_budget(configuration, requested):
    maximum = configuration["solver_wall_seconds_per_phase"]
    if requested is None:
        return maximum
    require(
        type(requested) in (int, float), "Solver wall budget must be numeric, not boolean/string"
    )
    try:
        value = float(requested)
    except OverflowError as error:
        raise ValueError("Solver wall budget is not finite") from error
    require(
        math.isfinite(value) and 60.0 <= value <= maximum,
        "Solver wall budget outside declared60..maximum range",
    )
    return value


def build_problem(
    candidate, ephemeris_path, *, phase, mission_profile=None, solver_wall_budget_seconds=None
):
    """Build one smooth graph. No solver construction/call or propagation."""
    _runtime()
    require(phase in PHASES, "Unknown local phase")
    profile = _identity(candidate, mission_profile)
    domain_contract = ephemeris.validate_domain(candidate.get("continuation_domain"))
    control_policy = controls.validate_policy(candidate.get("control_policy"))
    configuration = protocol(profile, domain_contract, control_policy)
    solver_budget = _solver_budget(configuration, solver_wall_budget_seconds)
    input_phase = "baseline" if phase == "restore" else "quality"
    duration, nodes, commands, target = validation.checked_candidate(candidate, input_phase)
    if phase == "quality":
        source_phase, source_status = candidate.get("solver_phase"), candidate.get("solver_status")
        valid_source = (
            (source_phase == "restore" and source_status == "Solve_Succeeded")
            or (
                source_phase == "quality"
                and source_status in {"Solve_Succeeded", *refinement.BUDGET_STATUSES}
            )
            or (source_phase == "handover" and source_status == "NOT_RUN_EXPLICIT_CELL_HANDOVER")
        )
        require(
            valid_source,
            "Quality phase requires a restoration, prior quality result or explicit handover",
        )
    selected = ephemeris.domain(candidate["arrival_cell_MJD"], domain_contract)
    arrival = build_arrival(ephemeris_path, selected["cell_MJD"], domain_contract)
    require(
        np.array_equal(target, np.asarray(arrival["evaluate"](duration)).reshape(6)),
        "Candidate target differs from exact arrival graph",
    )
    problem = model.build_problem(
        nodes[0],
        nodes,
        commands,
        transfer.LIGHTNESS,
        arrival,
        _settings(selected["cell_MJD"], domain_contract),
        duration,
    )
    w = problem["nlp"]["x"]
    physical = problem["functions"]["endpoint_physical"](w)
    restoration_scale = ca.DM(
        [profile["restoration_position_scale_m"] / model.AU_METRES] * 3
        + [profile["restoration_velocity_scale_m_s"] * model.DAY_SECONDS / model.AU_METRES] * 3
    )
    endpoint = physical / restoration_scale
    restoration_cost = 0.5 * ca.dot(endpoint, endpoint)
    # First restore a usable 600-day flight; only the next phase tries to reduce time.
    objective_expression = (
        restoration_cost
        if phase == "restore"
        else w[-1] * HORIZON / objective.dynamical_time_days(profile)
    )
    # Keep numerical room so rounding/replay does not push an on-limit answer outside the cap.
    optimisation_scale = ca.DM(
        [(profile["position_target_m"] - profile["position_nlp_reserve_m"]) / model.AU_METRES] * 3
        + [
            (profile["velocity_target_m_s"] - profile["velocity_nlp_reserve_m_s"])
            * model.DAY_SECONDS
            / model.AU_METRES
        ]
        * 3
    )
    qualified_endpoint = physical / optimisation_scale
    acceptance_scale = ca.DM(
        [profile["position_target_m"] / model.AU_METRES] * 3
        + [profile["velocity_target_m_s"] * model.DAY_SECONDS / model.AU_METRES] * 3
    )
    local_acceptance_endpoint = physical / acceptance_scale
    encounter = ca.vertcat(
        ca.dot(qualified_endpoint[:3], qualified_endpoint[:3]),
        ca.dot(qualified_endpoint[3:], qualified_endpoint[3:]),
    )
    constraints = ca.vertcat(problem["nlp"]["g"][: 9 * HOLDS], encounter)
    if phase == "restore":
        constraints = constraints[: 9 * HOLDS]
        problem["lbg"], problem["ubg"] = problem["lbg"][: 9 * HOLDS], problem["ubg"][: 9 * HOLDS]
        problem["lbx"][-1] = problem["ubx"][-1] = 1.0
        problem["layout"] = dict(problem["layout"], constraints=9 * HOLDS, free_variables=8 * HOLDS)
        problem["layout"].pop("encounter_norms_squared")
        problem["solver_options"]["ipopt.fixed_variable_treatment"] = "make_parameter"
        problem["limitations"].update(
            encounter_is_hard_constraint=False,
            fixed_time_feasibility_restoration_only=True,
            window_restricted_local_time_minimisation=False,
        )
    else:
        problem["limitations"].update(
            window_restricted_local_time_minimisation=True, minimum_time_at_declared_quality=True
        )
    guard = controls.symbolic_constraints(
        ca.reshape(w[6 * HOLDS : 8 * HOLDS], 2, HOLDS),
        w[-1] * HORIZON,
        HOLDS,
        control_policy,
        maximum_duration_days=HORIZON,
    )
    if control_policy is not None:
        start = int(constraints.numel())
        require(guard.numel() == HOLDS - 1, "Wrong adjacent-command guard count")
        constraints = ca.vertcat(constraints, guard)
        problem["lbg"] = np.r_[problem["lbg"], np.zeros(HOLDS - 1)]
        problem["ubg"] = np.r_[problem["ubg"], np.full(HOLDS - 1, np.inf)]
        problem["layout"] = dict(
            problem["layout"],
            constraints=start + HOLDS - 1,
            control_normal_step=[start, start + HOLDS - 1],
        )
    problem["control_regularity"] = {
        "policy": control_policy,
        "hard_constraints_applied": control_policy is not None,
        "physical_slew_certified": False,
    }
    problem["nlp"] = {"x": w, "f": objective_expression, "g": constraints}
    problem["functions"].update(
        f=ca.Function("fresh_" + phase + "_f", [w], [objective_expression]),
        g=ca.Function("fresh_" + phase + "_g", [w], [constraints]),
        gradient_f=ca.Function(
            "fresh_" + phase + "_gradient", [w], [ca.gradient(objective_expression, w)]
        ),
        jac_g=ca.Function("fresh_" + phase + "_jacobian", [w], [ca.jacobian(constraints, w)]),
        control_normal_step_margin=ca.Function("fresh_control_normal_step_margin", [w], [guard]),
    )
    global_acceptance_function = problem["functions"]["endpoint_acceptance_scaled"]
    problem["functions"].update(
        endpoint_scaled=ca.Function("fresh_local_endpoint", [w], [qualified_endpoint]),
        endpoint_norms_squared=ca.Function("fresh_local_encounter", [w], [encounter]),
        endpoint_acceptance_scaled=ca.Function(
            "fresh_local_acceptance_endpoint", [w], [local_acceptance_endpoint]
        ),
        endpoint_global_acceptance_scaled=global_acceptance_function,
        endpoint_restoration_scaled=ca.Function(
            "fresh_restoration_conditioning_endpoint", [w], [endpoint]
        ),
    )
    problem["numerical_reserve"] = {
        "position_m": profile["position_nlp_reserve_m"],
        "velocity_m_s": profile["velocity_nlp_reserve_m_s"],
        "NLP_position_limit_m": profile["position_target_m"] - profile["position_nlp_reserve_m"],
        "NLP_velocity_limit_m_s": profile["velocity_target_m_s"]
        - profile["velocity_nlp_reserve_m_s"],
        "local_position_target_m": profile["position_target_m"],
        "local_velocity_target_m_s": profile["velocity_target_m_s"],
        "hard_endpoint_constraints_applied": phase == "quality",
        "global_enc_limits_unchanged": True,
    }
    problem["solver_options"].update(
        {
            "ipopt.tol": 1e-10,
            "ipopt.constr_viol_tol": 1e-12,
            "ipopt.max_iter": configuration["iterations_per_phase"],
            "ipopt.max_wall_time": solver_budget,
            "ipopt.max_cpu_time": solver_budget,
            "ipopt.option_file_name": "",
        }
    )
    problem["phase"] = phase
    problem["arrival"] = arrival
    problem["protocol"] = configuration
    problem["solver_wall_budget_seconds"] = solver_budget
    return problem


def derivative_check(problem):
    """Bounded directional differences, including objective; no target extrapolation."""
    original = problem["w0"].copy()
    lower, upper = problem["lbx"], problem["ubx"]
    fixed = lower == upper
    probe = original.copy()
    probe[~fixed] = np.clip(probe[~fixed], lower[~fixed] + 2e-5, upper[~fixed] - 2e-5)
    rng = np.random.default_rng(600)
    checks = []
    for direction_index in range(3):
        direction = rng.normal(size=probe.size)
        direction[fixed] = 0.0
        direction /= np.linalg.norm(direction)
        analytic_g = np.asarray(problem["functions"]["jac_g"](probe) @ direction).reshape(-1)
        analytic_f = float(
            np.asarray(problem["functions"]["gradient_f"](probe)).reshape(-1) @ direction
        )
        for h in (1e-5, 3e-6, 1e-6):
            plus, minus = probe + h * direction, probe - h * direction
            require(
                np.all(plus >= lower)
                and np.all(plus <= upper)
                and np.all(minus >= lower)
                and np.all(minus <= upper),
                "Derivative stencil left bounds",
            )
            approximate_g = (
                np.asarray(problem["functions"]["g"](plus)).reshape(-1)
                - np.asarray(problem["functions"]["g"](minus)).reshape(-1)
            ) / (2 * h)
            approximate_f = (
                float(problem["functions"]["f"](plus)) - float(problem["functions"]["f"](minus))
            ) / (2 * h)
            error_g = np.abs(approximate_g - analytic_g)
            error_f = abs(approximate_f - analytic_f)
            passed = bool(
                np.isfinite(approximate_g).all()
                and math.isfinite(approximate_f)
                and np.all(error_g <= 2e-7 + 2e-6 * np.abs(analytic_g))
                and error_f <= 2e-7 + 2e-6 * abs(analytic_f)
            )
            checks.append(
                {
                    "direction": direction_index,
                    "step": h,
                    "constraints_max_error": float(error_g.max()),
                    "objective_error": error_f,
                    "passed": passed,
                }
            )
    require(np.array_equal(original, problem["w0"]), "Derivative check changed starting vector")
    return {
        "checks": checks,
        "passed": all(item["passed"] for item in checks),
        "fixed_variables_excluded": int(fixed.sum()),
        "seed_unchanged": True,
        "bounded_probes": True,
        "objective_checked": True,
        "scope": "sampled directional checks, not all-state derivative proof",
    }


def solve_candidate(
    candidate,
    ephemeris_path,
    *,
    phase,
    before_solver_call,
    mission_profile=None,
    solver_wall_budget_seconds=None,
):
    """One call and continuous host review; budget feasibility is not convergence.

    Parent must save the whole bundle and the unique marker made by the hook.
    The caller admits each same-run predecessor; this function rechecks its
    continuous physical history before each call. A parent
    still owns the outer worker deadline and the final native replay.
    """
    require(callable(before_solver_call), "Durable before_solver_call hook required")
    function_started = time.monotonic()
    print("Local step: Building the trajectory optimisation problem.", flush=True)
    seed = copy.deepcopy(candidate)
    started = time.monotonic()
    problem = build_problem(
        seed,
        ephemeris_path,
        phase=phase,
        mission_profile=mission_profile,
        solver_wall_budget_seconds=solver_wall_budget_seconds,
    )
    graph_seconds = time.monotonic() - started
    profile = problem["protocol"]["mission_profile"]
    # Recheck each physical start, including explicitly retimed previous-cell handovers.
    if phase == "quality":
        started = time.monotonic()
        starting_history_review = validation.validate_candidate(
            seed, ephemeris_path, stage="quality", mission_profile=profile
        )
        starting_history_recheck_seconds = time.monotonic() - started
        require(
            starting_history_review["accepted_for_native_review"],
            "Starting history is not independently host-accepted",
        )
    else:
        starting_history_review = None
        starting_history_recheck_seconds = None
    started = time.monotonic()
    print("Local step: Checking objective and constraint derivatives.", flush=True)
    derivatives = derivative_check(problem)
    require(derivatives["passed"], "Derivative check failed; solver not permitted")
    derivative_seconds = time.monotonic() - started
    started = time.monotonic()
    solver = ca.nlpsol("fresh_" + phase, "ipopt", problem["nlp"], problem["solver_options"])
    construction_seconds = time.monotonic() - started
    attempt = {
        "schema": SCHEMA,
        "phase": phase,
        "maximum_solver_calls": 1,
        "run_id": seed["run_id"],
        "origin_controller_sha256": seed["origin_controller_sha256"],
        "protocol": problem["protocol"],
        "options": dict(problem["solver_options"]),
        "solver_wall_budget_seconds": problem["solver_wall_budget_seconds"],
        "solver_iteration_limit": problem["solver_options"]["ipopt.max_iter"],
        "arrival_cell_MJD": list(problem["arrival"]["cell"]),
        "duration_bounds_days": list(problem["arrival"]["duration_bounds_days"]),
    }
    before_solver_call(copy.deepcopy(attempt))
    arguments = {key: problem[key] for key in ("lbx", "ubx", "lbg", "ubg")}
    started_wall, started_cpu = time.monotonic(), time.process_time()
    print("Local step: Solving with Ipopt. Trial values are not accepted trajectories.", flush=True)
    answer = solver(x0=problem["w0"], **arguments)
    wall_seconds, cpu_seconds = time.monotonic() - started_wall, time.process_time() - started_cpu
    stats = solver.stats()
    status, success, iterations = (
        stats.get(key) for key in ("return_status", "success", "iter_count")
    )
    metadata = (
        type(status) is str
        and type(success) is bool
        and type(iterations) is int
        and 0 <= iterations <= problem["solver_options"]["ipopt.max_iter"]
    )
    converged = metadata and status == "Solve_Succeeded" and success is True
    budget_eligible = bool(
        metadata
        and phase == "quality"
        and success is False
        and status in refinement.BUDGET_STATUSES
    )
    vector = np.asarray(answer["x"], dtype=float).reshape(-1)
    finite = vector.shape == (8 * HOLDS + 1,) and bool(np.isfinite(vector).all())
    bounds = finite and bool(np.all(vector >= problem["lbx"]) and np.all(vector <= problem["ubx"]))
    result_candidate, host, violation, returned_objective = None, None, None, None
    host_validation_seconds = None
    initial_objective = float(problem["functions"]["f"](problem["w0"]))
    errors = []
    if finite:
        nodes, commands, duration = model.unpack_solution(
            vector, STATE_SCALE, holds=HOLDS, duration_scale_days=HORIZON
        )
        result_candidate = dict(
            seed,
            nodes_cartesian=np.vstack([validation.INITIAL[:6], nodes]).tolist(),
            controls_clock_cone_rad=commands.tolist(),
            terminal_command_clock_cone_rad=commands[-1].tolist(),
            duration_days=duration,
            solver_status=status,
            solver_phase=phase,
            solver_target_cartesian=None,
            arrival_cell_MJD=list(problem["arrival"]["cell"]),
            duration_bounds_days=list(problem["arrival"]["duration_bounds_days"]),
        )
        try:
            lower_time, upper_time = problem["arrival"]["duration_bounds_days"]
            require(lower_time <= duration <= upper_time, "Returned duration outside cell")
            result_candidate["solver_target_cartesian"] = problem["arrival"]["evaluate"](
                duration
            ).tolist()
            g = np.asarray(problem["functions"]["g"](vector)).reshape(-1)
            require(
                np.isfinite(g).all() and g.shape == problem["lbg"].shape,
                "Invalid returned constraints",
            )
            returned_objective = float(problem["functions"]["f"](vector))
            require(math.isfinite(returned_objective), "Invalid returned objective")
            violation = max(
                0.0, float(np.max(problem["lbg"] - g)), float(np.max(g - problem["ubg"]))
            )
            # Check a budget-returned point too; physical feasibility and KKT convergence differ.
            if (
                (converged or budget_eligible)
                and bounds
                and violation <= 1e-12
                and returned_objective <= initial_objective + 1e-10
            ):
                started = time.monotonic()
                print(
                    "Local step: Checking the returned schedule by continuous propagation.",
                    flush=True,
                )
                host = validation.validate_candidate(
                    result_candidate, ephemeris_path, stage=phase, mission_profile=profile
                )
                host_validation_seconds = time.monotonic() - started
        except Exception as error:
            errors.append(type(error).__name__ + ": " + str(error))
    objective_not_increased = (
        returned_objective is not None and returned_objective <= initial_objective + 1e-10
    )
    eligible = bool(
        (converged or budget_eligible)
        and not errors
        and bounds
        and violation is not None
        and violation <= 1e-12
        and objective_not_increased
        and host is not None
        and host.get("accepted_for_native_review") is True
    )
    bundle = {
        "candidate": result_candidate,
        "result": {
            "schema": SCHEMA,
            "phase": phase,
            "status": (
                "HOST_CANDIDATE_NATIVE_PENDING" if eligible else "NO_ACCEPTED_LOCAL_CANDIDATE"
            ),
            "eligible_for_native_review": eligible,
            "solver_calls": 1,
            "solver_constructions": 1,
            "solver_wall_budget_seconds": problem["solver_wall_budget_seconds"],
            "solver_iteration_limit": problem["solver_options"]["ipopt.max_iter"],
            "solver_status": status if type(status) is str else None,
            "solver_reported_success": success if type(success) is bool else None,
            "solver_iterations": iterations if type(iterations) is int else None,
            "complete_convergence_metadata": bool(metadata),
            "solver_converged": bool(converged),
            "admission_policy": refinement.ADMISSION_POLICY,
            "admission_class": (
                ("converged" if converged else "budget_feasible_not_converged")
                if eligible
                else "rejected"
            ),
            "local_optimality_established": False,
            "finite_returned_vector": bool(finite),
            "original_bounds_met": bool(bounds),
            "constraint_violation": violation,
            "candidate_retained": result_candidate is not None,
            "initial_objective": initial_objective,
            "returned_objective": returned_objective,
            "objective_not_increased": objective_not_increased,
            "objective_comparison_tolerance": 1e-10,
            "native_mission_accepted": False,
            "new_accepted_improvement": False,
            "global_optimum_claim": False,
            "errors": errors,
        },
        "options": dict(problem["solver_options"]),
        "attempt": attempt,
        "derivative_check": derivatives,
        "host_validation": host,
        "starting_history_recheck": starting_history_review,
        "starting_history_phase": seed.get("solver_phase"),
        "timing": {
            "graph_seconds": graph_seconds,
            "derivative_seconds": derivative_seconds,
            "solver_construction_seconds": construction_seconds,
            "solver_wall_seconds": wall_seconds,
            "solver_cpu_seconds": cpu_seconds,
            "starting_history_recheck_seconds": starting_history_recheck_seconds,
            "host_validation_seconds": host_validation_seconds,
            "function_wall_seconds": time.monotonic() - function_started,
            "clock_relationship": (
                "function wall includes the reported substage clocks; do not sum twice"
            ),
        },
    }
    json.dumps(bundle, allow_nan=False)
    return bundle
