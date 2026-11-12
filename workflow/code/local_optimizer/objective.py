"""Arrival requirements, the local time objective and the native endpoint score.

These are different quantities. A score alone does not accept a trajectory.
"""

import copy
import math
import sys

SCHEMA = "mission_objective_v2"
DEFINITION = "validated_quality_then_minimum_time_v2"
AU_METRES = 149597870691.0
DAY_SECONDS = 86400.0
MU_SUN_M3_S2 = 1.32712440018e20
INITIAL_POSITION_AU = (-0.41386069061018826, 0.8923032332626103, -7.822802607859726e-06)
INITIAL_RADIUS_AU = math.sqrt(sum(value * value for value in INITIAL_POSITION_AU))
GATES = ("native_validated", "numerical_validated", "physical_validated")
REFERENCE_576_POSITION_M = 62990003.63941674
REFERENCE_576_VELOCITY_M_S = 63.027193526496035
BENCHMARK_STATUS = "measured_native_reference_quality_not_capture_requirement"


def default_profile():
    return {
        "schema": SCHEMA,
        "definition": DEFINITION,
        "mode": "minimum_time_at_declared_quality",
        "name": "earth_mercury_local_quality_development_v2",
        "position_target_m": 1000.0,
        "velocity_target_m_s": 0.01,
        "position_nlp_reserve_m": 250.0,
        "velocity_nlp_reserve_m_s": 0.0025,
        "restoration_position_scale_m": 100000000.0,
        "restoration_velocity_scale_m_s": 100.0,
        "global_position_limit_m": 100000000.0,
        "global_velocity_limit_m_s": 100.0,
        "numerical_position_gate_m": 100.0,
        "numerical_velocity_gate_m_s": 0.001,
        "fine_mesh_position_gate_m": 1.0,
        "fine_mesh_velocity_gate_m_s": 0.00001,
        "reference_radius_m": INITIAL_RADIUS_AU * AU_METRES,
        "solar_mu_m3_s2": MU_SUN_M3_S2,
        "target_status": "engineering_development_goal_not_capture_requirement",
    }


def reference_benchmark_profile():
    """No worse endpoint quality than the frozen 576-day native reference.

    These measured caps define a comparison, not a safe planetary approach.
    Keep the old development profile available for reproducing v2 results.
    """
    p = default_profile()
    p.update(
        name="earth_mercury_reference_quality_v1",
        position_target_m=REFERENCE_576_POSITION_M,
        velocity_target_m_s=REFERENCE_576_VELOCITY_M_S,
        target_status=BENCHMARK_STATUS,
    )
    return p


def _mission_number(value, name, *, positive=True):
    if type(value) not in (int, float):
        raise ValueError(name + " must be numeric, not boolean/string")
    try:
        value = float(value)
    except OverflowError as error:
        raise ValueError(name + " must be finite and representable") from error
    if not math.isfinite(value) or value < 0.0 or (positive and value == 0.0):
        raise ValueError(
            name + " must be finite and positive"
            if positive
            else name + " must be finite and nonnegative"
        )
    return value


def _checked_profile(profile):
    """Keep the checked profile and its time scale together for this call."""
    expected = default_profile()
    if profile is None:
        return expected, _time_scale(expected["reference_radius_m"], expected["solar_mu_m3_s2"])
    if type(profile) is not dict or set(profile) != set(expected):
        raise ValueError("Mission profile must contain exactly all v2 fields")
    result = copy.deepcopy(profile)
    for name in ("schema", "definition", "mode"):
        if result[name] != expected[name]:
            raise ValueError("Unsupported mission profile " + name)
    if result["target_status"] not in (expected["target_status"], BENCHMARK_STATUS):
        raise ValueError("Unsupported mission target status")
    if result["target_status"] == BENCHMARK_STATUS and any(
        result[k] != v
        for k, v in (
            ("position_target_m", REFERENCE_576_POSITION_M),
            ("velocity_target_m_s", REFERENCE_576_VELOCITY_M_S),
        )
    ):
        raise ValueError("Reference-quality caps must match the frozen native reference")
    if type(result["name"]) is not str or not result["name"].strip():
        raise ValueError("Mission profile name must be nonempty")
    for name in set(expected) - {"schema", "definition", "mode", "target_status", "name"}:
        result[name] = _mission_number(result[name], name)
    for prefix in ("position", "velocity"):
        unit = "m" if prefix == "position" else "m_s"
        target = result[prefix + "_target_" + unit]
        reserve = result[prefix + "_nlp_reserve_" + unit]
        numerical = result["numerical_" + prefix + "_gate_" + unit]
        fine = result["fine_mesh_" + prefix + "_gate_" + unit]
        if (
            not 2 * numerical + fine
            <= reserve
            < target
            <= result["global_" + prefix + "_limit_" + unit]
        ):
            raise ValueError(
                "Require two numerical gates plus fine mesh "
                "<= reserve < local target <= global limit"
            )
    # This version changes mission quality, never relaxes validated integration.
    if (
        result["numerical_position_gate_m"] != 100.0
        or result["numerical_velocity_gate_m_s"] != 0.001
    ):
        raise ValueError("Reviewed numerical gates are fixed at100m/.001m/s")
    if result["fine_mesh_position_gate_m"] != 1.0 or result["fine_mesh_velocity_gate_m_s"] != 1e-5:
        raise ValueError("Reviewed fine-mesh gates are fixed at1m/.00001m/s")
    scale = _time_scale(result["reference_radius_m"], result["solar_mu_m3_s2"])
    return result, scale


def validate_profile(profile=None):
    """Reject missing fields and unversioned changes; goals/scales are explicit."""
    if profile is None:
        return default_profile()
    return _checked_profile(profile)[0]


def _time_scale(radius: float, mu: float) -> float:
    # Radius and mu are checked at the caller's input boundary.
    # This order avoids overflowing radius cubed before the division.
    value = radius * math.sqrt(radius / mu) / DAY_SECONDS
    if not math.isfinite(value) or value <= 0.0:
        raise ValueError("Physical time scale is not finite and positive")
    return value


def dynamical_time_days(profile=None):
    """sqrt(r0^3/mu), in days; not the orbital period (no2pi factor)."""
    p = default_profile() if profile is None else profile
    radius = _mission_number(p.get("reference_radius_m"), "reference radius")
    mu = _mission_number(p.get("solar_mu_m3_s2"), "solar mu")
    return _time_scale(radius, mu)


def time_component(tof_days, profile=None):
    """Unvalidated preference only; it never establishes trajectory acceptance."""
    _, scale = _checked_profile(profile)
    duration = _mission_number(tof_days, "flight time", positive=False)
    return _time_record(duration, scale)


def _time_record(duration: float, scale: float):
    cost = duration / scale
    if not math.isfinite(cost):
        raise ValueError("Dimensionless time is not finite and representable")
    return {
        "definition": DEFINITION,
        "time_scale_days": scale,
        "time_scale_formula": "sqrt(reference_radius_m^3/solar_mu_m3_s2)/86400",
        "dimensionless_time_cost": cost,
        "unvalidated_time_score": 100.0 / (1.0 + cost),
        "validated_fitness": None,
        "acceptance_established": False,
    }


def evaluate(
    position_miss_m,
    velocity_miss_m_s,
    tof_days,
    profile=None,
    *,
    native_validated=None,
    numerical_validated=None,
    physical_validated=None,
):
    """Feasibility first; unknown validation stays null, never optimistic.

    Gate values must come from the caller's bound evidence. This pure function
    cannot authenticate a native receipt or prove a continuous trajectory.
    """
    p, scale = _checked_profile(profile)
    position = _mission_number(position_miss_m, "position miss", positive=False)
    velocity = _mission_number(velocity_miss_m_s, "velocity miss", positive=False)
    values = dict(
        native_validated=native_validated,
        numerical_validated=numerical_validated,
        physical_validated=physical_validated,
    )
    if any(value is not None and type(value) is not bool for value in values.values()):
        raise ValueError("Validation gates must be true, false or unknown")
    duration = _mission_number(tof_days, "flight time", positive=False)
    time = _time_record(duration, scale)
    quality = position <= p["position_target_m"] and velocity <= p["velocity_target_m_s"]
    known = all(value is not None for value in values.values())
    eligible = known and all(values.values()) and quality
    score = None if not known else time["unvalidated_time_score"] if eligible else 0.0
    benchmark = p["target_status"] == BENCHMARK_STATUS
    return {
        "schema": SCHEMA,
        "definition": DEFINITION,
        "profile": p,
        "position_miss_m": position,
        "velocity_miss_m_s": velocity,
        "tof_days": duration,
        "local_quality_limits_met": quality,
        "global_arrival_limits_met": position <= p["global_position_limit_m"]
        and velocity <= p["global_velocity_limit_m_s"],
        "validation_known": known,
        "fully_validated_quality_eligible": bool(eligible),
        "reference_benchmark_duration_days": 576.0 if benchmark else None,
        "reference_benchmark_time_improved": (
            bool(eligible and duration < 576.0 - 1e-6) if benchmark and known else None
        ),
        "gates": values,
        "own_fitness": score,
        "higher_is_better": True,
        **time,
        "validated_fitness": score,
        "acceptance_established": False,
        "scope": (
            "time ordering among quality-eligible validated trajectories; "
            "no global optimality claim"
        ),
    }


# Keep the native score contract distinct from the local mission profile.
NATIVE_CONTRACT_SCHEMA = "trajectory_score_contract_v1"
NATIVE_OWN_DEFINITION = "normalized_time_endpoint_squared_v1"
NATIVE_DEFINITION = "intrance_rendezvous_transfer_time_j_and_endpoint_v1"
NATIVE_EPS = 10.0 * sys.float_info.epsilon


def default_native_contract():
    """Return independent, explicit fixed scales; never normalize to each run."""
    return {
        "schema": NATIVE_CONTRACT_SCHEMA,
        "native_definition": NATIVE_DEFINITION,
        "own_definition": NATIVE_OWN_DEFINITION,
        "horizon_days": 600.0,
        "arrival_position_limit_m": 100_000_000.0,
        "arrival_velocity_limit_m_s": 100.0,
        "own_time_weight": 1.0,
        "own_position_weight": 0.5,
        "own_velocity_weight": 0.5,
    }


def _score_number(value, name, *, positive=False):
    if type(value) not in (int, float):
        raise ValueError(name + " must be a number, not a boolean or string")
    try:
        converted = float(value)
    except (ValueError, OverflowError) as exc:
        raise ValueError(name + " must be finite and representable") from exc
    if not math.isfinite(converted) or converted < 0 or (positive and converted == 0):
        raise ValueError(
            name
            + (" must be positive and finite" if positive else " must be nonnegative and finite")
        )
    return converted


def validate_native_contract(contract=None):
    """Reject incomplete/custom-unversioned contracts instead of guessing."""
    expected = default_native_contract()
    if contract is None:
        return expected
    if type(contract) is not dict or set(contract) != set(expected):
        raise ValueError("An explicit score contract must contain exactly all v1 fields")
    result = dict(contract)
    for name in ("schema", "native_definition", "own_definition"):
        if result[name] != expected[name]:
            raise ValueError("Unsupported " + name)
    # Positive weights keep each declared preference active; another objective
    # family needs a new definition, not a silent zero-weight switch.
    for name in set(expected) - {"schema", "native_definition", "own_definition"}:
        result[name] = _score_number(result[name], name, positive=True)
    return result


def _normalise_score(position_miss_m, velocity_miss_m_s, tof_days, contract):
    values = (
        _score_number(position_miss_m, "position_miss_m"),
        _score_number(velocity_miss_m_s, "velocity_miss_m_s"),
        _score_number(tof_days, "tof_days"),
    )
    ratios = (
        values[0] / contract["arrival_position_limit_m"],
        values[1] / contract["arrival_velocity_limit_m_s"],
        values[2] / contract["horizon_days"],
    )
    if not all(math.isfinite(value) for value in ratios):
        raise ValueError("Normalized metrics must be finite and representable")
    return values, ratios


def _native_components(rho_r, rho_v, normalized_time):
    # Keep the source's operation order; the simplified feasible formula hides
    # rounding at the boundary and the epsilon/large-error behavior.
    pos = math.log10(1.0 / max(rho_r, NATIVE_EPS))
    vel = math.log10(1.0 / max(rho_v, NATIVE_EPS))
    time = 100.0 * (1.0 - normalized_time)
    feasible = pos >= 0.0 and vel >= 0.0
    base = time if feasible else 0.0

    def recovered_error(subfitness):
        ratio = 1.0 / (1e-10 if subfitness < -10.0 else math.pow(10.0, subfitness))
        return 0.99 + ratio / 100.0 if ratio < 1.0 else ratio

    worst = max(recovered_error(pos), recovered_error(vel))
    score = 10.0 * base + 1.0 / math.sqrt(2.0 * worst * worst)
    if not all(math.isfinite(value) for value in (pos, vel, time, score)):
        raise ValueError("Native score is not finite and representable")
    return {
        "fitness": score,
        "position_subfitness": pos,
        "velocity_subfitness": vel,
        "time_subfitness": time,
        "time_reward_enabled": feasible,
    }


def native_endpoint_fitness(position_miss_m, velocity_miss_m_s, tof_days, *, contract=None):
    """Native-equivalent scalar at this endpoint, not a validation verdict."""
    contract = validate_native_contract(contract)
    _, ratios = _normalise_score(position_miss_m, velocity_miss_m_s, tof_days, contract)
    return _native_components(*ratios)["fitness"]
