"""Compare adjacent sail commands against the reference command changes.

Clock alpha and cone beta give n=[cos(beta), sin(beta)cos(alpha),
sin(beta)sin(alpha)] in the native RTN frame. The comparison uses one common
frame; it is a command jump, not actuator motion. The held commands still jump.
"""

import csv
import math
from pathlib import Path

import numpy as np

SCHEMA = "reference_control_regularity_v1"
# Measured from575consumed-command transitions in the precise576-day reference,
# using acos(clipped dot). The final metadata command is not another held arc.
REFERENCE_MAX_STEP_RAD = 0.16302811923339847
ROUNDING_RAD = 1e-10


def frozen_reference_policy():
    return {
        "schema": SCHEMA,
        "mode": "reference_relative_command_grid",
        "reference_id": "native_576d_precise_controller",
        "reference_time_step_days": 1.0,
        "max_step_rad": REFERENCE_MAX_STEP_RAD,
        "max_step_rate_rad_per_day": REFERENCE_MAX_STEP_RAD,
        "angular_roundoff_allowance_rad": ROUNDING_RAD,
        "physical_slew_certified": False,
    }


def validate_policy(policy=None):
    """Absent means the legacy unguarded contract; this v1 has one fixed reference."""
    if policy is None:
        return None
    expected = frozen_reference_policy()
    if type(policy) is not dict or set(policy) != set(expected):
        raise ValueError("Control policy must contain exactly all declared v1 fields")
    for key, value in expected.items():
        supplied = policy[key]
        if type(value) is float:
            if (
                type(supplied) not in (int, float)
                or not math.isfinite(supplied)
                or supplied != value
            ):
                raise ValueError("Control reference cap/roundoff changed: " + key)
        elif type(supplied) is not type(value) or supplied != value:
            raise ValueError("Unsupported control policy field: " + key)
    return dict(expected)


def _controls(controls):
    raw = np.asarray(controls, dtype=object)
    if (
        raw.ndim != 2
        or raw.shape[1] != 2
        or len(raw) < 2
        or not all(
            isinstance(v, (int, float, np.integer, np.floating))
            and not isinstance(v, (bool, np.bool_))
            for v in raw.flat
        )
    ):
        raise ValueError("Expected at least two held clock/cone pairs, not terminal metadata")
    values = np.asarray(controls, dtype=float)
    if (
        not np.isfinite(values).all()
        or np.any(np.abs(values[:, 0]) > math.pi)
        or np.any((values[:, 1] < 0.0) | (values[:, 1] > math.pi / 2.0))
    ):
        raise ValueError("Nonfinite or out-of-domain physical command")
    return values


def normal_vectors(controls):
    return _normal_vectors(_controls(controls))


def _normal_vectors(values: np.ndarray) -> np.ndarray:
    # The public caller has already checked these clock/cone pairs.
    clock, cone = values[:, 0], values[:, 1]
    return np.column_stack(
        (np.cos(cone), np.sin(cone) * np.cos(clock), np.sin(cone) * np.sin(clock))
    )


def diagnose(controls, duration_days, policy=None):
    """Inspect held commands on their declared uniform grid; never propagate.

    The rate is jump divided by the preceding hold duration, a grid proxy only.
    It is not a maximum instantaneous physical angular rate.
    """
    values = _controls(controls)
    if (
        type(duration_days) not in (int, float)
        or not math.isfinite(duration_days)
        or duration_days <= 0.0
    ):
        raise ValueError("Control duration must be finite and positive")
    p = validate_policy(policy)
    normals = _normal_vectors(values)
    dots = np.sum(normals[:-1] * normals[1:], axis=1)
    angles = np.arccos(np.clip(dots, -1.0, 1.0))
    hold = float(duration_days) / len(values)
    if not math.isfinite(hold) or hold <= 0.0:
        raise ValueError("Unrepresentable control hold duration")
    rate = angles / hold
    if not np.isfinite(rate).all():
        raise ValueError("Unrepresentable command-grid rate")
    cap = None if p is None else min(p["max_step_rad"], p["max_step_rate_rad_per_day"] * hold)
    passed = (
        None if p is None else bool(np.all(angles <= cap + p["angular_roundoff_allowance_rad"]))
    )
    return {
        "schema": SCHEMA,
        "enabled": p is not None,
        "policy": p,
        "passed": passed,
        "held_intervals": len(values),
        "command_transitions": len(values) - 1,
        "duration_days": float(duration_days),
        "hold_duration_days": hold,
        "max_step_rad": float(angles.max()),
        "max_step_deg": float(np.rad2deg(angles.max())),
        "max_step_rate_rad_per_day": float(rate.max()),
        "max_step_rate_deg_per_day": float(np.rad2deg(rate.max())),
        "total_command_variation_rad": float(angles.sum()),
        "allowed_step_rad": cap,
        "angular_roundoff_allowance_rad": (
            None if p is None else p["angular_roundoff_allowance_rad"]
        ),
        "maximum_step_excess_rad": None if cap is None else float(max(0.0, angles.max() - cap)),
        "physical_slew_certified": False,
        "inertial_frame_motion_included": False,
        "terminal_metadata_excluded": True,
        "zero_order_hold_discontinuities_removed": False,
    }


def build_policy(reference_commands_path):
    """Recompute the frozen cap from the exact reference command file.

    The caller pins the file identity in its private provenance record. This
    function checks the reference grid, consumed rows and measured cap; no
    absolute path or file fingerprint is embedded in the public policy.
    """
    path = Path(reference_commands_path)
    if not path.is_file() or path.is_symlink():
        raise ValueError("Use the exact regular reference command file")
    with path.open(newline="", encoding="ascii") as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames != ["MJD_d", "clock_rad", "cone_rad"]:
            raise ValueError("Reference command columns changed")
        rows = list(reader)
    if len(rows) != 577 or any(set(row) != {"MJD_d", "clock_rad", "cone_rad"} for row in rows):
        raise ValueError("Reference requires576holds and577boundaries")
    times = np.array([float(row["MJD_d"]) for row in rows])
    if not np.array_equal(times, 52654.5 + np.arange(577)):
        raise ValueError("Reference departure/one-day timing changed")
    commands = [[float(row["clock_rad"]), float(row["cone_rad"])] for row in rows]
    _controls(commands)  # Terminal metadata must remain physically well formed.
    report = diagnose(commands[:-1], 576.0)
    expected = frozen_reference_policy()
    # Only floating arithmetic roundoff is tolerated, never a different cap.
    if abs(report["max_step_rad"] - expected["max_step_rad"]) > ROUNDING_RAD:
        raise ValueError("Measured reference normal-step maximum changed")
    return expected


def symbolic_constraints(controls, total_days, holds, policy, *, maximum_duration_days=600.0):
    """Smooth dot-minus-cos residuals >=0, with no angular inverse in the NLP.

    The time-scaled cap is tighter than the absolute cap throughout the checked
    duration domain. This avoids a min() kink at the600-day upper boundary.
    """
    import casadi as ca

    p = validate_policy(policy)
    if p is None:
        return ca.MX.zeros(0, 1)
    if type(holds) is not int or holds < 2 or tuple(controls.shape) != (2, holds):
        raise ValueError("Symbolic controls must have shape2xN")
    if (
        type(maximum_duration_days) not in (int, float)
        or not math.isfinite(maximum_duration_days)
        or not (0 < maximum_duration_days <= holds * p["reference_time_step_days"])
    ):
        raise ValueError("Rate-only symbolic form would not enforce the absolute reference cap")
    clock, cone = controls[0, :], controls[1, :]
    # Same-frame normals make 179 to -179 degrees a small turn, not a 358-degree jump.
    dots = ca.cos(cone[:-1]) * ca.cos(cone[1:]) + ca.sin(cone[:-1]) * ca.sin(cone[1:]) * ca.cos(
        clock[1:] - clock[:-1]
    )
    threshold = p["max_step_rate_rad_per_day"] * total_days / holds
    return ca.vec(dots - ca.cos(threshold))
