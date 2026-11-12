"""Read physical native records; units and column order are part of the interface."""

import csv
import math

DEPARTURE_MJD = 52654.5


LIGHTNESS = 0.09274704077328881


COMMAND_FIELDS = ("MJD_d", "clock_rad", "cone_rad")


INITIAL_FIELDS = tuple(
    (
        "MJD_d,r_AU,phi_rad,theta_rad,dr_AU_per_d,dphi_rad_per_d,dtheta_rad_per_d,"
        "mp_kg,srd_native,x_AU,y_AU,z_AU,vx_AU_per_d,vy_AU_per_d,vz_AU_per_d,"
        "lightness,ideal_sail_flag"
    ).split(",")
)


STATE_FIELDS = (
    ("command_index",)
    + INITIAL_FIELDS[:15]
    + ("clock_rad", "cone_rad", "lightness", "ideal_sail_flag")
)


LIMIT_FIELDS = tuple(
    ("min_duration_d,max_duration_d,max_interval_d,min_radius_AU,max_radius_AU").split(",")
)


CARTESIAN_FIELDS = INITIAL_FIELDS[9:15]


FIXED_INITIAL = [
    DEPARTURE_MJD,
    0.9836085254736493,
    2.0050765264748667,
    -7.953166737930806e-06,
    5.558569973125507e-05,
    0.017771459999688165,
    -5.691039321939054e-07,
    0.0,
    0.0,
    -0.41386069061018826,
    0.8923032332626103,
    -7.822802607859726e-06,
    -0.015880919317195288,
    -0.007304482858971418,
    -5.602175619070018e-07,
    LIGHTNESS,
    1.0,
]


MODEL = {
    "propulsion": "ideal_solar_sail",
    "gravity": "solar_only",
    "state_frame": "heliocentric_native_J2000_ecliptic",
    "control_frame": "instantaneous_radial_transverse_normal",
    "control_convention": "n=cos(cone)*er+sin(cone)*(cos(clock)*et+sin(clock)*eh)",
    "target": "mercury",
    "ephemeris": "DE405",
    "native_profile": "mercury_body_rendezvous_v1",
}


UNITS = {
    "epoch": "MJD_day",
    "position": "AU",
    "velocity": "AU/day",
    "control": "rad",
    "duration": "day",
}


ENDPOINT_FIELDS = tuple(
    (
        "requested_final_MJD_d,actual_final_MJD_d,duration_d,reached_end,"
        "completed_intervals,command_rows,min_sampled_radius_AU,max_sampled_radius_AU,"
        "max_boundary_clock_error_d,path_checks_passed,initial_state_verified,"
        "x_AU,y_AU,z_AU,vx_AU_per_d,vy_AU_per_d,vz_AU_per_d,final_mp_kg,"
        "final_srd_native,diagnostic_only"
    ).split(",")
)


TARGET_FIELDS = tuple(
    (
        "requested_final_MJD_d,actual_final_MJD_d,target_body_id,target_MJD_d,"
        "spacecraft_x_AU,spacecraft_y_AU,spacecraft_z_AU,spacecraft_vx_AU_per_d,"
        "spacecraft_vy_AU_per_d,spacecraft_vz_AU_per_d,target_x_AU,target_y_AU,"
        "target_z_AU,target_vx_AU_per_d,target_vy_AU_per_d,target_vz_AU_per_d,"
        "position_residual_AU,velocity_residual_AU_per_d,diagnostic_only"
    ).split(",")
)


INTERVAL_FIELDS = tuple(
    (
        "interval_index,start_MJD_d,requested_end_MJD_d,actual_end_MJD_d,duration_d,"
        "RK_status,reached_boundary,mp_kg,srd_native,radius_AU"
    ).split(",")
)


def _require(condition, message):
    if not condition:
        raise ValueError(message)


def _read_numeric(path, fields, limit=1000):
    with path.open(newline="", encoding="ascii") as stream:
        reader = csv.reader(stream)
        _require(next(reader, None) == list(fields), "CSV fields/units mismatch: " + path.name)
        rows = []
        for line_number, row in enumerate(reader, 2):
            _require(
                len(rows) < limit and len(row) == len(fields),
                "CSV count/width mismatch: " + path.name,
            )
            try:
                values = [float(value) for value in row]
            except ValueError as error:
                raise ValueError(f"Invalid CSV number at {path.name}:{line_number}") from error
            _require(
                all(math.isfinite(value) for value in values), "nonfinite CSV value: " + path.name
            )
            rows.append(dict(zip(fields, values)))
    _require(rows, "empty CSV: " + path.name)
    return rows


def _angles(clock, cone):
    _require(
        -math.pi <= clock <= math.pi and 0.0 <= cone <= math.pi / 2,
        "clock/cone outside native angle bounds",
    )


def _check_polar(row):
    radius, phi, theta, dr, dphi, dtheta = [row[k] for k in INITIAL_FIELDS[1:7]]
    _require(
        0.2 <= radius <= 1.5 and abs(math.cos(theta)) > 1e-8,
        "sampled state outside radius/latitude domain",
    )
    ct, st, cp, sp = math.cos(theta), math.sin(theta), math.cos(phi), math.sin(phi)
    radial = [ct * cp, ct * sp, st]
    azimuthal = [-sp, cp, 0.0]
    latitudinal = [-st * cp, -st * sp, ct]
    predicted = [radius * x for x in radial] + [
        dr * radial[j] + radius * ct * dphi * azimuthal[j] + radius * dtheta * latitudinal[j]
        for j in range(3)
    ]
    # This checks representation roundoff, not integration or mission accuracy.
    _require(
        all(
            abs(a - row[k]) <= 2e-14 * max(1.0, abs(a)) for a, k in zip(predicted, CARTESIAN_FIELDS)
        ),
        "polar/Cartesian representation mismatch",
    )
