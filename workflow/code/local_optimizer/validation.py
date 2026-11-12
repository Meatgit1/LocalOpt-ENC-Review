"""Check one returned Mercury candidate and its native replay.

Arrival miss, shooting continuity and numerical agreement are separate checks.
Continuous propagation starts from the original state, without node resets.
"""

from pathlib import Path
import csv
import math
import re

import numpy as np

from . import controls, ephemeris, objective, propagation, records, transfer
from .transfer import ENDPOINT_FIELDS, INTERVAL_FIELDS, TARGET_FIELDS

HOLDS = 600
CELL = (53248.0, 53256.0)
DURATION = (593.5, 600.0)
MESHES = (64, 128, 256)
PROFILE = "mercury_body_rendezvous_precision_v3"
SCALE = np.array([1.0, 1.0, 1.0, 0.02, 0.02, 0.02])
INITIAL = np.array(transfer.FIXED_INITIAL[9:15] + transfer.FIXED_INITIAL[7:9])
JET_LIMITS = (1e-3, 1e-7, 1e-10, 1e-13)  # m, m/s, m/s², m/s³


def require(ok, message):
    if not ok:
        raise ValueError(message)


def _number(value, name):
    require(type(value) in (int, float) and math.isfinite(value), "Invalid " + name)
    return float(value)


def _array(value, shape, name):
    raw = np.asarray(value, dtype=object)
    require(
        raw.shape == shape
        and all(
            isinstance(x, (int, float, np.integer, np.floating))
            and not isinstance(x, (bool, np.bool_))
            for x in raw.flat
        ),
        "Invalid shape/type: " + name,
    )
    result = np.asarray(value, dtype=float)
    require(np.isfinite(result).all(), "Nonfinite " + name)
    return result


def _stage(stage):
    require(
        stage in ("baseline", "restore", "restoration", "quality", "handover"),
        "Unknown fresh600 stage",
    )
    return "restore" if stage == "restoration" else stage


def _profile(candidate=None, supplied=None):
    source = (
        supplied
        if supplied is not None
        else candidate.get("mission_profile") if candidate is not None else None
    )
    profile = objective.validate_profile(source)
    original = objective.default_profile()
    for key in (
        "reference_radius_m",
        "solar_mu_m3_s2",
        "global_position_limit_m",
        "global_velocity_limit_m_s",
    ):
        require(profile[key] == original[key], "Mercury physics/global setting changed: " + key)
    if candidate is not None:
        require(
            candidate.get("mission_profile") == profile,
            "Candidate/validation mission profile differs",
        )
    return profile


def contract(mission_profile=None, *, cell_MJD=None, continuation_domain=None, control_policy=None):
    profile = _profile(supplied=mission_profile)
    selected = ephemeris.domain(cell_MJD, continuation_domain=continuation_domain)
    return {
        "schema": "fresh600_validation_v4",
        "mission_profile": profile,
        "holds": HOLDS,
        "state_count": HOLDS + 1,
        "departure_MJD": transfer.DEPARTURE_MJD,
        "arrival_cell_MJD": selected["cell_MJD"],
        "duration_bounds_days": selected["duration_bounds_days"],
        "independent_meshes": list(MESHES),
        "cell_domain_contract": ephemeris.SCHEMA,
        "continuation_domain": ephemeris.validate_domain(continuation_domain),
        "control_policy": controls.validate_policy(control_policy),
        "local_continuity_RK_substeps": 16,
        "continuity_scaled_limit": 1e-9,
        "fine_mesh_position_m": 1.0,
        "fine_mesh_velocity_m_s": 1e-5,
        "native_agreement_position_m": 100.0,
        "native_agreement_velocity_m_s": 0.001,
        "position_limit_km": profile["position_target_m"] / 1000.0,
        "velocity_limit_m_s": profile["velocity_target_m_s"],
        "global_position_limit_km": profile["global_position_limit_m"] / 1000.0,
        "global_velocity_limit_m_s": profile["global_velocity_limit_m_s"],
        "radius_bounds_AU": [0.2, 1.5],
        "velocity_component_limit_AU_day": 0.1,
        "minimum_angular_momentum_AU2_day": 1e-8,
        "minimum_cos_latitude": 1e-8,
        "minimum_transverse_rate_per_day": 1e-12,
        "target_jet_limits_SI": list(JET_LIMITS),
        "minimum_time_reduction_days": 1e-6,
        "target_epoch_shift_max_seconds": 2.05e-5,
        "target_epoch_alignment_position_m": 0.001,
        "target_epoch_alignment_velocity_m_s": 1e-7,
        "target_epoch_linear_residual_position_m": 0.002,
        "sampled_path_only": True,
        "independent_ephemeris_data_source": False,
    }


def checked_candidate(candidate, stage="baseline", *, mission_profile=None):
    stage = _stage(stage)
    require(
        type(candidate) is dict
        and type(candidate.get("schema")) is int
        and candidate["schema"] == 1,
        "Candidate schema",
    )
    _profile(candidate, mission_profile)
    controls.validate_policy(candidate.get("control_policy"))
    require(
        type(candidate.get("holds")) is int
        and candidate["holds"] == HOLDS
        and type(candidate.get("state_count")) is int
        and candidate["state_count"] == HOLDS + 1,
        "Candidate fresh600 layout must be600holds/601states",
    )
    require(
        _number(candidate.get("departure_MJD"), "departure") == transfer.DEPARTURE_MJD
        and _number(candidate.get("lightness"), "lightness") == transfer.LIGHTNESS,
        "Original departure/model changed",
    )
    duration_bounds = ephemeris.bounds(candidate)
    duration = _number(candidate.get("duration_days"), "duration")
    require(
        duration_bounds[0] <= duration <= duration_bounds[1]
        and (
            stage in ("quality", "handover")
            or (tuple(candidate["arrival_cell_MJD"]) == CELL and duration == 600.0)
        ),
        "Duration outside stage/cell contract",
    )
    nodes = _array(candidate.get("nodes_cartesian"), (HOLDS + 1, 6), "nodes")
    control_angles = _array(candidate.get("controls_clock_cone_rad"), (HOLDS, 2), "controls")
    target = _array(candidate.get("solver_target_cartesian"), (6,), "solver target")
    require(np.array_equal(nodes[0], INITIAL[:6]), "Original initial state changed")
    require(
        np.all(np.abs(control_angles[:, 0]) <= math.pi)
        and np.all((control_angles[:, 1] >= 0) & (control_angles[:, 1] <= math.pi / 2)),
        "Sail control bounds",
    )
    terminal = _array(candidate.get("terminal_command_clock_cone_rad"), (2,), "terminal command")
    transfer._angles(float(terminal[0]), float(terminal[1]))
    require(
        stage == "baseline" or np.array_equal(terminal, control_angles[-1]),
        "Restoration/time terminal metadata must repeat the last held command",
    )
    return duration, nodes, control_angles, target


def _terminal_command(candidate, supplied):
    terminal = np.asarray(candidate["terminal_command_clock_cone_rad"], dtype=float)
    if supplied is not None:
        require(
            np.array_equal(terminal, _array(supplied, (2,), "terminal command")),
            "Explicit terminal differs from candidate",
        )
    return terminal


def _command_rows(values, terminal):
    duration, _, controls_array, _ = values
    times = [transfer.DEPARTURE_MJD + duration * i / HOLDS for i in range(HOLDS + 1)]
    times[-1] = transfer.DEPARTURE_MJD + duration
    return [
        dict(MJD_d=t, clock_rad=float(u[0]), cone_rad=float(u[1]))
        for t, u in zip(times, np.vstack((controls_array, terminal)))
    ]


def expected_commands(candidate, *, stage="baseline", terminal_command=None, mission_profile=None):
    values = checked_candidate(candidate, stage, mission_profile=mission_profile)
    return _command_rows(values, _terminal_command(candidate, terminal_command))


def write_schedule(candidate, output, *, stage="baseline", mission_profile=None):
    """Write exactly one fresh schedule directory; preserve existing attempts."""
    commands = expected_commands(candidate, stage=stage, mission_profile=mission_profile)
    output = Path(output)
    require(not output.exists(), "Schedule directory exists; no overwrite")
    attained = commands[-1]["MJD_d"] - transfer.DEPARTURE_MJD
    rows = {
        "bridge_initial_input.csv": (transfer.INITIAL_FIELDS, [transfer.FIXED_INITIAL]),
        "bridge_input.csv": (
            transfer.COMMAND_FIELDS,
            [[row[key] for key in transfer.COMMAND_FIELDS] for row in commands],
        ),
        "bridge_schedule_limits.csv": (
            transfer.LIMIT_FIELDS,
            [[attained, attained, 2.0, 0.2, 1.5]],
        ),
    }
    output.mkdir(parents=True)
    for name, (fields, values) in rows.items():
        with (output / name).open("x", newline="", encoding="ascii") as stream:
            writer = csv.writer(stream, lineterminator="\n")
            writer.writerow(fields)
            writer.writerows([[format(float(value), ".17g") for value in row] for row in values])
    _check_schedule_rows(output, commands, candidate["duration_days"])
    return output


def check_schedule(output, candidate, *, stage="baseline", mission_profile=None):
    """Read the schedule again at its file boundary."""
    expected = expected_commands(candidate, stage=stage, mission_profile=mission_profile)
    return _check_schedule_rows(Path(output), expected, candidate["duration_days"])


def _check_schedule_rows(output, expected, requested_duration):
    required = ("bridge_input.csv", "bridge_initial_input.csv", "bridge_schedule_limits.csv")
    require(
        all((output / name).is_file() and not (output / name).is_symlink() for name in required),
        "Missing/linked schedule file",
    )
    commands = transfer._read_numeric(output / required[0], transfer.COMMAND_FIELDS, HOLDS + 1)
    initial = transfer._read_numeric(output / required[1], transfer.INITIAL_FIELDS, 1)
    limits = transfer._read_numeric(output / required[2], transfer.LIMIT_FIELDS, 1)
    attained = expected[-1]["MJD_d"] - transfer.DEPARTURE_MJD
    require(
        commands == expected
        and len(initial) == len(limits) == 1
        and list(initial[0].values()) == transfer.FIXED_INITIAL
        and list(limits[0].values()) == [attained, attained, 2.0, 0.2, 1.5],
        "Schedule controls, timing, original initial state or limits differ",
    )
    return {
        "schedule_checked": True,
        "holds": HOLDS,
        "duration_days": attained,
        "requested_duration_days": requested_duration,
        "native_acceptance": False,
    }


def _error(first, second):
    a, b = np.asarray(first), np.asarray(second)
    require(
        a.shape == b.shape
        and a.ndim == 2
        and a.shape[1] >= 6
        and np.isfinite(a).all()
        and np.isfinite(b).all(),
        "Invalid comparison arrays",
    )
    delta = a - b
    return {
        "position_m": float(np.max(np.linalg.norm(delta[:, :3], axis=1))) * propagation.AU_METRES,
        "velocity_m_s": float(np.max(np.linalg.norm(delta[:, 3:6], axis=1)))
        * propagation.AU_METRES
        / propagation.DAY_SECONDS,
    }


def independent_target(duration, ephemeris_path, cell_MJD=None, continuation_domain=None):
    """Direct NumPy differentiated polynomials, plus separately rounded Clenshaw.

    Coefficients are selected by the cell's interior midpoint, not by copying the
    CasADi graph. Raw data are shared; arithmetic/derivative paths are independent.
    """
    selected = ephemeris.domain(cell_MJD, continuation_domain=continuation_domain)
    cell = selected["cell_MJD"]
    duration_bounds = selected["duration_bounds_days"]
    duration = _number(duration, "target duration")
    require(duration_bounds[0] <= duration <= duration_bounds[1], "Target extrapolation refused")
    eph = ephemeris.Ephemeris(ephemeris_path)
    record = eph.records[eph.record_index(sum(cell) / 2 + ephemeris.JD0)]
    pieces = []
    acceleration_bound = 0.0
    for sign, body in ((1.0, 0), (-1.0, 10)):
        entry, count, groups = eph.pointers[body]
        width = float((record[1] - record[0]) / groups)
        group = int(((sum(cell) / 2 + ephemeris.JD0) - record[0]) / width)
        start = float(record[0] - ephemeris.JD0 + group * width)
        require(
            0 <= group < groups and start <= cell[0] < cell[1] <= start + width,
            "Raw Mercury/Sun interval does not contain the declared cell",
        )
        offset = entry - 1 + 3 * count * group
        coeff = record[offset : offset + 3 * count].reshape(3, count)
        pieces.append((sign, start, width, coeff))
        axis_bound = np.array(
            [sum(np.abs(np.polynomial.chebyshev.chebder(axis, m=2))) for axis in coeff]
        )
        acceleration_bound += (
            float(np.linalg.norm(axis_bound)) * (2.0 / width / ephemeris.DAY_S) ** 2 * 1000.0
        )

    def orders(value, *, rounded_native_epoch=False):
        # Native JD rounding may lie a few microseconds outside a clipped NLP
        # window. It still must remain inside this same raw polynomial cell.
        allowed = (
            selected["coefficient_duration_bounds_days"]
            if rounded_native_epoch
            else duration_bounds
        )
        require(
            allowed[0] <= value <= allowed[1], "Rounded target outside chosen raw coefficient cell"
        )
        result = [np.zeros(3) for _ in range(4)]
        for sign, start, width, coeff in pieces:
            # Small offsets express departure+T-start. Use the raw cell start,
            # never the lower bound of a clipped admissible duration window.
            cell_start_duration = cell[0] - transfer.DEPARTURE_MJD
            tau = 2.0 * ((cell[0] - start) + (value - cell_start_duration)) / width - 1.0
            require(-1.0 <= tau <= 1.0, "Polynomial extrapolation")
            for order in range(4):
                result[order] += (
                    sign
                    * np.array(
                        [
                            np.polynomial.chebyshev.chebval(
                                tau, np.polynomial.chebyshev.chebder(axis, m=order)
                            )
                            for axis in coeff
                        ]
                    )
                    * (2.0 / width) ** order
                )
        return [ephemeris.ROT @ q / (ephemeris.AU_M / 1000.0) for q in result]

    jet = orders(duration)
    smooth = np.r_[jet[0], jet[1]]
    rounded = eph.heliocentric("mercury", transfer.DEPARTURE_MJD + duration)
    native_jd = (transfer.DEPARTURE_MJD + duration) + ephemeris.JD0
    native_duration = (native_jd - ephemeris.JD0) - transfer.DEPARTURE_MJD
    aligned_jet = orders(native_duration, rounded_native_epoch=True)
    aligned = np.r_[aligned_jet[0], aligned_jet[1]]
    shift = (native_duration - duration) * ephemeris.DAY_S
    alignment = _error(aligned[None, :], rounded[None, :])
    linear_residual = float(
        np.linalg.norm(
            (rounded[:3] - smooth[:3]) * ephemeris.AU_M
            - smooth[3:] * ephemeris.AU_M / ephemeris.DAY_S * shift
        )
    )
    require(
        all(np.isfinite(q).all() for q in jet) and np.isfinite(rounded).all(), "Nonfinite target"
    )
    return {
        "state": smooth,
        "derivative": np.r_[jet[1], jet[2]],
        "second_derivative": np.r_[jet[2], jet[3]],
        "native_epoch_state": rounded,
        "arrival_cell_MJD": list(cell),
        "duration_bounds_days": list(duration_bounds),
        "aligned_native_epoch": alignment,
        "epoch_shift_seconds": shift,
        "linear_shift_residual_m": linear_residual,
        "acceleration_bound_m_s2": acceleration_bound,
    }


def audit_target_join(ephemeris_path, first_cell, second_cell, continuation_domain=None):
    """Audit both one-sided polynomials at a reviewed join and nearby rounded epochs.

    Position/velocity discrepancies must fit existing target arithmetic gates.
    Higher derivative jumps are reported, not silently assumed to vanish.
    This is an ephemeris check, not a spacecraft or native-process replay.
    """
    pair = sorted(
        [
            ephemeris.domain(first_cell, continuation_domain=continuation_domain),
            ephemeris.domain(second_cell, continuation_domain=continuation_domain),
        ],
        key=lambda item: item["cell_MJD"][0],
    )
    left, right = pair
    require(
        left["cell_MJD"][1] == right["cell_MJD"][0], "Target cells are not distinct adjacent cells"
    )
    duration = left["cell_MJD"][1] - transfer.DEPARTURE_MJD
    left_value = independent_target(duration, ephemeris_path, left["cell_MJD"], continuation_domain)
    right_value = independent_target(
        duration, ephemeris_path, right["cell_MJD"], continuation_domain
    )
    jump = _error(left_value["state"][None, :], right_value["state"][None, :])
    require(
        jump["position_m"] <= 0.001 and jump["velocity_m_s"] <= 1e-7,
        "Target position/velocity join exceeds declared arithmetic gates",
    )
    probes = []
    for selected, sign in ((left, -math.inf), (right, math.inf)):
        for probe in (duration, math.nextafter(duration, sign)):
            value = independent_target(
                probe, ephemeris_path, selected["cell_MJD"], continuation_domain
            )
            aligned = value["aligned_native_epoch"]
            require(
                aligned["position_m"] <= 0.001
                and aligned["velocity_m_s"] <= 1e-7
                and abs(value["epoch_shift_seconds"]) <= 2.05e-5
                and value["linear_shift_residual_m"]
                <= 0.002
                + 0.5 * value["acceleration_bound_m_s2"] * value["epoch_shift_seconds"] ** 2,
                "Native-rounded target join is inconsistent",
            )
            probes.append(
                {
                    "cell_MJD": selected["cell_MJD"],
                    "duration_days": probe,
                    "epoch_shift_seconds": value["epoch_shift_seconds"],
                    "native_epoch_alignment": aligned,
                }
            )
    acceleration = (
        float(np.linalg.norm(left_value["derivative"][3:] - right_value["derivative"][3:]))
        * ephemeris.AU_M
        / ephemeris.DAY_S**2
    )
    jerk = (
        float(
            np.linalg.norm(
                left_value["second_derivative"][3:] - right_value["second_derivative"][3:]
            )
        )
        * ephemeris.AU_M
        / ephemeris.DAY_S**3
    )
    return {
        "schema": "fresh600_target_join_v1",
        "join_duration_days": duration,
        "left_cell_MJD": left["cell_MJD"],
        "right_cell_MJD": right["cell_MJD"],
        "state_jump": jump,
        "acceleration_jump_m_s2": acceleration,
        "jerk_jump_m_s3": jerk,
        "position_velocity_gate_passed": True,
        "one_sided_derivative_jumps_reported": True,
        "rounded_epoch_probes": probes,
        "global_smoothness_claimed": False,
        "native_process_executed": False,
        "spacecraft_propagation_performed": False,
    }


def validate_target_model(ephemeris_path, arrival=None, cell_MJD=None, continuation_domain=None):
    """Six bounded target/derivative probes; no optimisation or propagation."""
    if cell_MJD is None and arrival is not None:
        cell_MJD = arrival.get("cell")
    selected = ephemeris.domain(cell_MJD, continuation_domain=continuation_domain)
    cell, duration_bounds = selected["cell_MJD"], selected["duration_bounds_days"]
    sources = [
        Path(ephemeris_path),
        Path(__file__),
        Path(ephemeris.__file__),
    ]
    bindings = {str(path.resolve()): records.sha(path) for path in sources}
    descriptor = ephemeris.descriptor(ephemeris_path, cell, continuation_domain=continuation_domain)
    if arrival is None:
        arrival = ephemeris.build_cell(ephemeris_path, descriptor)
    require(
        tuple(arrival["cell"]) == tuple(cell)
        and tuple(arrival["duration_bounds_days"]) == tuple(duration_bounds),
        "Target graph cell/bounds differ",
    )
    probes = []
    lo, hi = duration_bounds
    width, offset = hi - lo, min(1e-6, (hi - lo) / 8.0)
    for duration in (lo, lo + offset, lo + 0.25 * width, (lo + hi) / 2.0, hi - offset, hi):
        direct = independent_target(duration, ephemeris_path, cell, continuation_domain)
        state = _array(arrival["evaluate"](duration), (6,), "graph target")
        first = _array(arrival["evaluate_derivative"](duration), (6,), "graph derivative")
        second = _array(
            arrival["evaluate_second_derivative"](duration), (6,), "graph second derivative"
        )
        # Check both appearances of v and a; all derivatives are with respect to DAYS.
        errors = [
            float(np.linalg.norm(state[:3] - direct["state"][:3])) * ephemeris.AU_M,
            max(
                float(np.linalg.norm(state[3:] - direct["state"][3:])),
                float(np.linalg.norm(first[:3] - direct["derivative"][:3])),
            )
            * ephemeris.AU_M
            / ephemeris.DAY_S,
            max(
                float(np.linalg.norm(first[3:] - direct["derivative"][3:])),
                float(np.linalg.norm(second[:3] - direct["second_derivative"][:3])),
            )
            * ephemeris.AU_M
            / ephemeris.DAY_S**2,
            float(np.linalg.norm(second[3:] - direct["second_derivative"][3:]))
            * ephemeris.AU_M
            / ephemeris.DAY_S**3,
        ]
        require(
            all(error <= limit for error, limit in zip(errors, JET_LIMITS)),
            "Target graph/independent derivative disagreement",
        )
        require(
            direct["aligned_native_epoch"]["position_m"] <= 0.001
            and direct["aligned_native_epoch"]["velocity_m_s"] <= 1e-7
            and abs(direct["epoch_shift_seconds"]) <= 2.05e-5
            and direct["linear_shift_residual_m"]
            <= 0.002 + 0.5 * direct["acceleration_bound_m_s2"] * direct["epoch_shift_seconds"] ** 2,
            "Target epoch alignment inconsistent",
        )
        probes.append(
            {
                "duration_days": duration,
                "jet_errors_SI": errors,
                "native_epoch_alignment": direct["aligned_native_epoch"],
            }
        )
    require(
        all(records.sha(path) == digest for path, digest in bindings.items()),
        "Target evidence changed during validation",
    )
    joins = [
        audit_target_join(ephemeris_path, cell, neighbour, continuation_domain)
        for neighbour in ephemeris.neighbours(cell, continuation_domain=continuation_domain)
    ]
    require(
        all(records.sha(path) == digest for path, digest in bindings.items()),
        "Target evidence changed during join validation",
    )
    return {
        "schema": "fresh600_target_cell_validation_v3",
        "target_cell_checked": True,
        "arrival_cell_MJD": list(cell),
        "duration_bounds_days": list(duration_bounds),
        "probes": probes,
        "continuation_domain": ephemeris.validate_domain(continuation_domain),
        "exact_cell_descriptor": descriptor,
        "join_audits": joins,
        "jet_limits_SI": list(JET_LIMITS),
        "direct_native_target_check_pending": True,
        "global_smoothness_claimed": False,
        "shared_raw_ephemeris": True,
        "evidence_sha256": bindings,
    }


def _domain(states, stats, mesh):
    require(
        states.shape == (HOLDS + 1, 8) and np.isfinite(states).all(),
        "Independent trajectory shape/finite failure",
    )
    require(
        stats["initializations"] == 1
        and stats["shooting_resets"] == 0
        and stats["completed_intervals"] == HOLDS
        and stats["rhs_calls"] == 4 * mesh * HOLDS
        and stats["substeps_per_interval"] == mesh,
        "Independent integration accounting",
    )
    require(
        0.2 <= stats["minimum_radius_AU"] <= stats["maximum_radius_AU"] <= 1.5
        and stats["minimum_angular_momentum_AU2_day"] >= 1e-8
        and stats["minimum_cosine_latitude"] > 1e-8
        and stats["minimum_transverse_angular_motion_per_day"] > 1e-12
        and stats["above_cutoff_stage_count"] == 0
        and stats["minimum_mass_kg"] == stats["maximum_mass_kg"] == 0.0
        and stats["minimum_exposure_native"] >= 0.0
        and np.all(np.diff(states[:, 7]) >= 0.0)
        and np.all(np.abs(states[:, 3:6]) <= 0.1),
        "Independent sampled path outside original domain",
    )
    require(np.array_equal(states[0], INITIAL), "Independent initial state changed")


def _continuous(values, terminal_command):
    commands = _command_rows(values, terminal_command)
    histories, statistics = {}, {}
    for mesh in MESHES:
        # Fly from departure again; never reset to the solver's intermediate nodes here.
        values, stats = propagation.propagate(
            INITIAL.tolist(), commands, propagation.Model(transfer.LIGHTNESS), mesh
        )
        values = np.asarray(values)
        _domain(values, stats, mesh)
        histories[mesh], statistics[str(mesh)] = values, stats
    return histories, statistics, commands


def _nodes_domain(nodes):
    radius = np.linalg.norm(nodes[:, :3], axis=1)
    momentum = np.linalg.norm(np.cross(nodes[:, :3], nodes[:, 3:6]), axis=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        cosine = np.linalg.norm(nodes[:, :2], axis=1) / radius
        transverse = momentum / radius**2
    return bool(
        np.all((radius >= 0.2) & (radius <= 1.5))
        and np.all(np.abs(nodes[:, 3:6]) <= 0.1)
        and np.all(momentum >= 1e-8)
        and np.all(cosine > 1e-8)
        and np.all(transverse > 1e-12)
    )


def _mission(state, target, mission_profile=None):
    profile = _profile(supplied=mission_profile)
    position = float(np.linalg.norm(state[:3] - target[:3])) * ephemeris.AU_M / 1000.0
    velocity = float(np.linalg.norm(state[3:6] - target[3:])) * ephemeris.AU_M / ephemeris.DAY_S
    return {
        "position_miss_km": position,
        "velocity_miss_m_s": velocity,
        "limits_pass": position * 1000.0 <= profile["position_target_m"]
        and velocity <= profile["velocity_target_m_s"],
        "global_limits_pass": position * 1000.0 <= profile["global_position_limit_m"]
        and velocity <= profile["global_velocity_limit_m_s"],
        "position_target_m": profile["position_target_m"],
        "velocity_target_m_s": profile["velocity_target_m_s"],
    }


def _local_checks(values, stage, histories):
    duration, nodes, control_angles, _ = values
    node_error = _error(nodes, histories[256][:, :6])
    defects = []
    if stage != "baseline":
        for index in range(HOLDS):
            state = nodes[index].tolist() + INITIAL[6:].tolist()
            for _ in range(16):
                state = propagation.rk4(
                    state,
                    duration / HOLDS / 16,
                    *control_angles[index],
                    propagation.Model(transfer.LIGHTNESS),
                )
            defects.append(float(np.max(np.abs((nodes[index + 1] - state[:6]) / SCALE))))
    node_domain_pass = _nodes_domain(nodes)
    passed = stage == "baseline" or (
        node_domain_pass
        and node_error["position_m"] <= 100.0
        and node_error["velocity_m_s"] <= 0.001
        and max(defects) <= 1e-9
    )
    return {
        "local_consistency_pass": bool(passed),
        "nodes_vs_continuous256": node_error,
        "maximum_local_RK16_continuity_scaled": max(defects) if defects else None,
        "node_continuity_gate_applied": stage != "baseline",
        "node_domain_pass": node_domain_pass,
    }


def validate_candidate(
    candidate, ephemeris_path, *, stage="baseline", terminal_command=None, mission_profile=None
):
    """Host gate: one initialisation per mesh; baseline infeasibility is allowed."""
    stage = _stage(stage)
    profile = _profile(candidate, mission_profile)
    values = checked_candidate(candidate, stage, mission_profile=profile)
    duration, _, control_angles, proposed_target = values
    cell = candidate["arrival_cell_MJD"]
    context = candidate.get("continuation_domain")
    target = independent_target(duration, ephemeris_path, cell, context)
    regularity = controls.diagnose(control_angles, duration, candidate.get("control_policy"))
    regularity_pass = not regularity["enabled"] or regularity["passed"] is True
    histories, statistics, _ = _continuous(values, _terminal_command(candidate, terminal_command))
    fine = _error(histories[128], histories[256])
    local = _local_checks(values, stage, histories)
    target_error = _error(proposed_target[None, :], target["state"][None, :])
    encounter = _mission(histories[256][-1], target["native_epoch_state"], profile)
    epoch_consistent = (
        target["aligned_native_epoch"]["position_m"] <= 0.001
        and target["aligned_native_epoch"]["velocity_m_s"] <= 1e-7
        and abs(target["epoch_shift_seconds"]) <= 2.05e-5
        and target["linear_shift_residual_m"]
        <= 0.002 + 0.5 * target["acceleration_bound_m_s2"] * target["epoch_shift_seconds"] ** 2
    )
    numeric = (
        epoch_consistent
        and fine["position_m"] <= 1.0
        and fine["velocity_m_s"] <= 1e-5
        and target_error["position_m"] <= 0.001
        and target_error["velocity_m_s"] <= 1e-7
        and local["local_consistency_pass"]
    )
    return {
        "schema": "fresh600_candidate_validation_v4",
        "stage": stage,
        "validation_contract": contract(
            profile,
            cell_MJD=cell,
            continuation_domain=context,
            control_policy=candidate.get("control_policy"),
        ),
        "control_regularity": regularity,
        "control_regularity_pass": regularity_pass,
        "host_numerical_accepted": bool(numeric),
        "mission_limits_pass_without_numerical_gate": encounter["limits_pass"],
        "accepted_for_native_review": bool(
            numeric and (stage == "baseline" or (encounter["limits_pass"] and regularity_pass))
        ),
        "fine_mesh128_vs256": fine,
        **local,
        "target_comparison": target_error,
        "target_epoch_alignment_pass": bool(epoch_consistent),
        "independent_statistics": statistics,
        "encounter": encounter,
        "duration_days": duration,
        "native_numerical_accepted": False,
        "native_mission_accepted": False,
        "native_process_executed": False,
        "local_solve_performed_by_validator": False,
        "time_reduction_days": 600.0 - duration,
        "time_reduction_above_resolution": duration < 600.0 - 1e-6,
        "time_reduction_verified": False,
        "lineage_verified": False,
        "sampled_path_only": True,
        "global_optimum_claimed": False,
    }


def _read(directory, name, fields, count):
    rows = transfer._read_numeric(Path(directory) / name, fields, count)
    require(len(rows) == count, "Native record count: " + name)
    return rows


def inspect_native_records(
    fixture,
    result,
    candidate,
    ephemeris_path,
    *,
    stage="baseline",
    terminal_command=None,
    mission_profile=None,
):
    """Scientific file binding only; process-receipt acceptance belongs to transport."""
    fixture, result = Path(fixture), Path(result)
    stage = _stage(stage)
    profile = _profile(candidate, mission_profile)
    values = checked_candidate(candidate, stage, mission_profile=profile)
    duration, _, _, proposed_target = values
    sim = (fixture / "mercury.sim").read_text()
    for key in ("MAX_RELATIVE_ERROR", "MAX_ABSOLUTE_ERROR"):
        require(
            re.findall(r"(?m)^" + key + r"\s*=\s*([^\s/]+)", sim) == ["1e-14"],
            "V3 precision must be exactly1e-14",
        )
    stdout = (result / "stdout.txt").read_text()
    require(
        re.findall(r"NATIVE_BRIDGE_PROFILE_COMPLETE profile=([^\s]+)", stdout) == [PROFILE]
        and stdout.count("NATIVE_SCHEDULE_REACHED_END") == 1,
        "Wrong/missing native V3 schedule completion",
    )
    expected = _command_rows(values, _terminal_command(candidate, terminal_command))
    commands = _read(fixture, "bridge_input.csv", transfer.COMMAND_FIELDS, HOLDS + 1)
    require(
        commands == expected
        and _read(result, "bridge_commands.csv", transfer.COMMAND_FIELDS, HOLDS + 1) == commands,
        "Native schedule differs from exact controls/timing/terminal metadata",
    )
    first = _read(fixture, "bridge_initial_input.csv", transfer.INITIAL_FIELDS, 1)[0]
    require(
        list(first.values()) == transfer.FIXED_INITIAL
        and _read(result, "bridge_initial.csv", transfer.INITIAL_FIELDS, 1)[0] == first,
        "Native original initial state changed",
    )
    epoch = commands[-1]["MJD_d"]
    attained = epoch - transfer.DEPARTURE_MJD
    require(
        abs(attained - duration) <= math.ulp(epoch),
        "Native duration beyond epoch-rounding allowance",
    )
    require(
        list(_read(fixture, "bridge_schedule_limits.csv", transfer.LIMIT_FIELDS, 1)[0].values())
        == [attained, attained, 2.0, 0.2, 1.5],
        "Native schedule bounds changed",
    )
    rows = _read(result, "bridge_states.csv", transfer.STATE_FIELDS, HOLDS + 1)
    prior = first["srd_native"]
    for index, (row, command) in enumerate(zip(rows, commands)):
        require(
            row["command_index"] == index
            and all(row[key] == command[key] for key in transfer.COMMAND_FIELDS)
            and row["lightness"] == transfer.LIGHTNESS
            and row["ideal_sail_flag"] == 1.0
            and row["mp_kg"] == 0.0
            and row["srd_native"] >= prior,
            "Native boundary/model accounting",
        )
        transfer._check_polar(row)
        prior = row["srd_native"]
    require(
        [rows[0][key] for key in transfer.INITIAL_FIELDS] == transfer.FIXED_INITIAL,
        "Native first boundary differs",
    )
    intervals = _read(result, "bridge_intervals.csv", INTERVAL_FIELDS, HOLDS)
    for index, (record, left, right, row) in enumerate(
        zip(intervals, commands, commands[1:], rows[1:])
    ):
        require(
            record["interval_index"] == index
            and record["start_MJD_d"] == left["MJD_d"]
            and record["requested_end_MJD_d"] == record["actual_end_MJD_d"] == right["MJD_d"]
            and record["duration_d"] == right["MJD_d"] - left["MJD_d"]
            and record["RK_status"] == record["reached_boundary"] == 1.0
            and record["mp_kg"] == row["mp_kg"]
            and record["srd_native"] == row["srd_native"]
            and record["radius_AU"] == row["r_AU"],
            "Native interval incomplete/inconsistent",
        )
    end = _read(result, "bridge_endpoint.csv", ENDPOINT_FIELDS, 1)[0]
    encounter = _read(result, "bridge_target_endpoint.csv", TARGET_FIELDS, 1)[0]
    require(
        end["completed_intervals"] == HOLDS
        and end["command_rows"] == HOLDS + 1
        and end["requested_final_MJD_d"] == end["actual_final_MJD_d"] == epoch
        and end["duration_d"] == attained
        and all(
            end[key] == 1.0
            for key in (
                "reached_end",
                "path_checks_passed",
                "initial_state_verified",
                "diagnostic_only",
            )
        ),
        "Native endpoint accounting",
    )
    require(
        all(end[key] == rows[-1][key] for key in transfer.CARTESIAN_FIELDS)
        and end["final_mp_kg"] == rows[-1]["mp_kg"]
        and end["final_srd_native"] == rows[-1]["srd_native"],
        "Native endpoint state mismatch",
    )
    require(
        0.2 <= end["min_sampled_radius_AU"] <= end["max_sampled_radius_AU"] <= 1.5
        and 0.0 <= end["max_boundary_clock_error_d"] <= 1e-9,
        "Native sampled path/clock bound",
    )
    require(
        encounter["requested_final_MJD_d"]
        == encounter["actual_final_MJD_d"]
        == encounter["target_MJD_d"]
        == epoch
        and encounter["target_body_id"] == 2.0
        and encounter["diagnostic_only"] == 1.0
        and all(
            encounter["spacecraft_" + key] == rows[-1][key] for key in transfer.CARTESIAN_FIELDS
        ),
        "Native target epoch/body/spacecraft mismatch",
    )
    values = np.array(
        [
            [row[key] for key in transfer.CARTESIAN_FIELDS] + [row["mp_kg"], row["srd_native"]]
            for row in rows
        ]
    )
    target = np.array([encounter["target_" + key] for key in transfer.CARTESIAN_FIELDS])
    cell = candidate["arrival_cell_MJD"]
    context = candidate.get("continuation_domain")
    direct = independent_target(attained, ephemeris_path, cell, context)
    smooth = independent_target(duration, ephemeris_path, cell, context)["state"]
    regularity = controls.diagnose(
        [[row["clock_rad"], row["cone_rad"]] for row in commands[:-1]],
        attained,
        candidate.get("control_policy"),
    )
    regularity_pass = not regularity["enabled"] or regularity["passed"] is True
    require(
        stage == "baseline" or regularity_pass,
        "Native held commands fail declared reference regularity",
    )
    solver_target_error = _error(proposed_target[None, :], smooth[None, :])
    require(
        solver_target_error["position_m"] <= 0.001 and solver_target_error["velocity_m_s"] <= 1e-7,
        "Candidate target differs from independently decoded smooth target",
    )
    alignment = _error(target[None, :], direct["native_epoch_state"][None, :])
    require(
        alignment["position_m"] <= 0.001 and alignment["velocity_m_s"] <= 1e-7,
        "Native target disagrees with same-epoch independentDE405",
    )
    dr, dv = np.linalg.norm(values[-1, :3] - target[:3]), np.linalg.norm(
        values[-1, 3:6] - target[3:]
    )
    require(
        math.isclose(float(dr), encounter["position_residual_AU"], rel_tol=1e-12, abs_tol=1e-18)
        and math.isclose(
            float(dv), encounter["velocity_residual_AU_per_d"], rel_tol=1e-12, abs_tol=1e-18
        ),
        "Native residuals disagree with vectors",
    )
    return {
        "candidate": candidate,
        "stage": stage,
        "terminal_command": terminal_command,
        "native": values,
        "native_target": target,
        "native_encounter": _mission(values[-1], target, profile),
        "mission_profile": profile,
        "control_regularity": regularity,
        "control_regularity_pass": regularity_pass,
        "native_cartesian_states": values[:, :6].copy(),
        "initial_MJD": transfer.DEPARTURE_MJD,
        "actual_final_MJD": epoch,
        "solver_target_comparison": solver_target_error,
        "target_alignment": alignment,
        "duration_days": attained,
        "commands": commands,
        "scientific_records_bound": True,
        "process_receipt_accepted": False,
        "profile": PROFILE,
        "tolerance": "1e-14",
        "initial": INITIAL.copy(),
    }


def numerical_review(info, output_directory=None, *, mission_profile=None):
    """Run the independently continuous mesh ladder, without guessed-node resets."""
    profile = _profile(info["candidate"], mission_profile)
    require(info.get("mission_profile") == profile, "Native review mission profile differs")
    if output_directory is not None:
        output_directory = Path(output_directory)
        require(not output_directory.exists(), "History output exists; preserve prior review")
    values = checked_candidate(info["candidate"], info["stage"], mission_profile=profile)
    histories, statistics, _ = _continuous(
        values, _terminal_command(info["candidate"], info["terminal_command"])
    )
    fine = _error(histories[128], histories[256])
    local = _local_checks(values, info["stage"], histories)
    differences = {str(mesh): _error(info["native"], histories[mesh]) for mesh in MESHES}
    numeric = bool(
        fine["position_m"] <= 1.0
        and fine["velocity_m_s"] <= 1e-5
        and all(
            value["position_m"] <= 100.0 and value["velocity_m_s"] <= 0.001
            for value in differences.values()
        )
        and local["local_consistency_pass"]
    )
    mission = info["native_encounter"]["limits_pass"]
    recorded_policy = info["candidate"].get("control_policy")
    if recorded_policy is not None:
        commands = info["commands"]
        repeated = controls.diagnose(
            [[row["clock_rad"], row["cone_rad"]] for row in commands[:-1]],
            info["duration_days"],
            recorded_policy,
        )
        require(
            info.get("control_regularity") == repeated
            and info.get("control_regularity_pass") is repeated["passed"],
            "Native command regularity verdict differs from bound command rows",
        )
    regularity_pass = info.get("control_regularity_pass", True)
    require(type(regularity_pass) is bool, "Malformed control regularity verdict")
    saved = {}
    if output_directory is not None:
        output_directory.mkdir(parents=True)
        for mesh, values in histories.items():
            path = output_directory / ("continuous_RK" + str(mesh) + ".npy")
            with path.open("xb") as stream:
                np.save(stream, values, allow_pickle=False)
            saved[str(mesh)] = {"path": str(path.resolve()), "sha256": records.sha(path)}
    return {
        "schema": "fresh600_native_numerical_v4",
        "stage": info["stage"],
        "native_numerical_accepted": False,
        "validation_contract": contract(
            profile,
            cell_MJD=info["candidate"]["arrival_cell_MJD"],
            continuation_domain=info["candidate"].get("continuation_domain"),
            control_policy=info["candidate"].get("control_policy"),
        ),
        "control_regularity": info.get("control_regularity"),
        "control_regularity_pass": regularity_pass,
        "recorded_native_numerical_checks_passed": numeric,
        "mission_limits_pass_without_numerical_gate": mission,
        "native_mission_accepted": False,
        "numerical_and_mission_gates_pass": bool(numeric and mission and regularity_pass),
        "duration_days": info["duration_days"],
        "time_reduction_days": 600.0 - info["duration_days"],
        "time_reduction_above_resolution": info["duration_days"] < 600.0 - 1e-6,
        "time_reduction_verified": False,
        "fine_mesh128_vs256": fine,
        "native_vs_continuous": differences,
        **local,
        "solver_target_comparison": info["solver_target_comparison"],
        "independent_statistics": statistics,
        "encounter": info["native_encounter"],
        "native_target_vs_same_epoch_DE405": info["target_alignment"],
        "scientific_records_bound": True,
        "process_receipt_accepted": False,
        "profile": PROFILE,
        "tolerance": "1e-14",
        "sampled_path_only": True,
        "saved_histories": saved,
        "no_shooting_resets": True,
        "global_optimum_claimed": False,
        "local_solve_performed_by_validator": False,
    }
