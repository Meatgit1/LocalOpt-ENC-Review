"""Independent Cartesian RK4 arithmetic for a bounded native sail bridge audit.

State: x,y,z [AU], vx,vy,vz [AU/day], native mass coordinate [kg],
native accumulated solar exposure. No ephemeris, optimiser or shooting resets.
"""

from dataclasses import dataclass
import math

AU_METRES = 149597870691.0
DAY_SECONDS = 86400.0
GM_SI = 1.32712440018e20
GM_AUD = GM_SI * DAY_SECONDS**2 / AU_METRES**3


def require(condition, reason):
    if not condition:
        raise ValueError(reason)


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def norm(a):
    return math.sqrt(dot(a, a))


def cross(a, b):
    return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]


@dataclass(frozen=True)
class Model:
    lightness: float
    mu: float = GM_AUD
    sail_cutoff_AU: float = 5.0

    def __post_init__(self):
        require(math.isfinite(self.lightness) and self.lightness >= 0, "Invalid exported lightness")
        require(math.isfinite(self.mu) and self.mu > 0, "Invalid gravitational constant")


def cartesian_from_polar(y):
    r, p, t, dr, dp, dt = y[:6]
    ct, st, cp, sp = math.cos(t), math.sin(t), math.cos(p), math.sin(p)
    er = [ct * cp, ct * sp, st]
    ep = [-sp, cp, 0.0]
    el = [-st * cp, -st * sp, ct]
    return (
        [r * q for q in er]
        + [dr * er[j] + r * ct * dp * ep[j] + r * dt * el[j] for j in range(3)]
        + list(y[6:8])
    )


def rhs(state, clock, cone, model, sampled=None):
    require(len(state) == 8 and all(math.isfinite(q) for q in state), "Nonfinite/invalid state")
    require(-math.pi <= clock <= math.pi and 0 <= cone <= math.pi / 2, "Sail angle bounds violated")
    r, v = state[:3], state[3:6]
    radius = norm(r)
    require(radius > 0, "Zero solar radius")
    h = cross(r, v)
    hn = norm(h)
    require(hn > 1e-15, "Degenerate RTN angular momentum")
    er = [q / radius for q in r]
    eh = [q / hn for q in h]
    et = cross(eh, er)
    ca, sa, cb, sb = math.cos(clock), math.sin(clock), math.cos(cone), math.sin(cone)
    sail_normal = [cb * er[j] + sb * (ca * et[j] + sa * eh[j]) for j in range(3)]
    factor = 0.0 if radius > model.sail_cutoff_AU else model.lightness * model.mu / radius**2
    force = [factor * cb * cb * q for q in sail_normal]
    acceleration = [-model.mu * r[j] / radius**3 + force[j] for j in range(3)]
    if sampled is not None:
        cosine_latitude = math.hypot(r[0], r[1]) / radius
        transverse_angular_motion = hn / radius**2
        sampled["rhs_calls"] += 1
        sampled["minimum_radius_AU"] = min(sampled["minimum_radius_AU"], radius)
        sampled["maximum_radius_AU"] = max(sampled["maximum_radius_AU"], radius)
        sampled["minimum_angular_momentum_AU2_day"] = min(
            sampled["minimum_angular_momentum_AU2_day"], hn
        )
        sampled["minimum_mass_kg"] = min(sampled["minimum_mass_kg"], state[6])
        sampled["maximum_mass_kg"] = max(sampled["maximum_mass_kg"], state[6])
        sampled["above_cutoff_stage_count"] += int(radius > model.sail_cutoff_AU)
        sampled["minimum_cosine_latitude"] = min(
            sampled["minimum_cosine_latitude"], cosine_latitude
        )
        sampled["minimum_transverse_angular_motion_per_day"] = min(
            sampled["minimum_transverse_angular_motion_per_day"], transverse_angular_motion
        )
        sampled["minimum_exposure_native"] = min(sampled["minimum_exposure_native"], state[7])
    # The native polar source disables thrust beyond 5 AU but not dose accumulation.
    value = list(v) + acceleration + [0.0, cb / radius**2]
    require(all(math.isfinite(q) for q in value), "Nonfinite independent RHS")
    return value


def rk4(state, h, clock, cone, model, sampled=None):
    def shifted(k, fraction):
        return [x + fraction * h * q for x, q in zip(state, k)]

    k1 = rhs(state, clock, cone, model, sampled)
    k2 = rhs(shifted(k1, 0.5), clock, cone, model, sampled)
    k3 = rhs(shifted(k2, 0.5), clock, cone, model, sampled)
    k4 = rhs(shifted(k3, 1.0), clock, cone, model, sampled)
    value = [x + h * (a + 2 * b + 2 * c + d) / 6 for x, a, b, c, d in zip(state, k1, k2, k3, k4)]
    require(all(math.isfinite(q) for q in value), "Nonfinite RK4 boundary state")
    return value


def propagate(initial, commands, model, substeps):
    """Initialize once; [time_i,time_(i+1)) holds use command_i.

    Last command is terminal metadata and has no propagation interval. Native
    boundary states are never inputs to this routine after initial construction.
    """
    require(isinstance(substeps, int) and substeps > 0, "Invalid fixed refinement")
    require(len(commands) >= 2, "Need at least one complete hold")
    for row in commands:
        require(set(row) == {"MJD_d", "clock_rad", "cone_rad"}, "Command schema mismatch")
        require(all(math.isfinite(q) for q in row.values()), "Nonfinite command")
        require(
            -math.pi <= row["clock_rad"] <= math.pi and 0 <= row["cone_rad"] <= math.pi / 2,
            "Control bound violation",
        )
    stats = {
        "rhs_calls": 0,
        "minimum_radius_AU": math.inf,
        "maximum_radius_AU": 0.0,
        "minimum_angular_momentum_AU2_day": math.inf,
        "minimum_mass_kg": math.inf,
        "maximum_mass_kg": -math.inf,
        "above_cutoff_stage_count": 0,
        "minimum_cosine_latitude": math.inf,
        "minimum_transverse_angular_motion_per_day": math.inf,
        "minimum_exposure_native": math.inf,
        "initializations": 1,
        "shooting_resets": 0,
        "completed_intervals": 0,
        "substeps_per_interval": substeps,
    }
    state = list(initial)
    states = [list(state)]
    for left, right in zip(commands, commands[1:]):
        duration = right["MJD_d"] - left["MJD_d"]
        require(duration > 0, "Commands not strictly increasing")
        step = duration / substeps
        for _ in range(substeps):
            state = rk4(state, step, left["clock_rad"], left["cone_rad"], model, stats)
        states.append(list(state))
        stats["completed_intervals"] += 1
    require(
        stats["rhs_calls"] == 4 * substeps * (len(commands) - 1),
        "Incomplete independent integration",
    )
    return states, stats


def differences(a, b):
    require(len(a) == len(b), "Boundary count mismatch")
    rows = []
    for i, (x, y) in enumerate(zip(a, b)):
        require(
            len(x) == len(y) == 8 and all(math.isfinite(q) for q in [*x, *y]),
            "Nonfinite/invalid comparison state",
        )
        rows.append(
            {
                "boundary": i,
                "position_m": norm([u - v for u, v in zip(x[:3], y[:3])]) * AU_METRES,
                "velocity_m_s": norm([u - v for u, v in zip(x[3:6], y[3:6])])
                * AU_METRES
                / DAY_SECONDS,
                "mass_kg": abs(x[6] - y[6]),
                "exposure_native": abs(x[7] - y[7]),
            }
        )
    return rows


def maxima(rows):
    return {
        key: max(row[key] for row in rows)
        for key in ("position_m", "velocity_m_s", "mass_kg", "exposure_native")
    }
