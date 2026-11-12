"""Mercury at the arrival time: numeric DE405 decoding and bounded cell graphs.

The symbolic graph uses one coefficient cell. The numeric decoder and NumPy jets
provide a separate arithmetic check; all use the same raw DE405 data.
"""

from __future__ import annotations

import copy
import math
from pathlib import Path
import struct

import casadi as ca
import numpy as np

RECORD_BYTES = 8144


JD0 = 2400000.5


AU_M = 149597870691.0


DAY_S = 86400.0


EPS_ARCSEC = 84381.412


EPS = EPS_ARCSEC * (math.pi / 180.0) / 3600.0


ROT = np.array(
    [[1.0, 0.0, 0.0], [0.0, math.cos(EPS), math.sin(EPS)], [0.0, -math.sin(EPS), math.cos(EPS)]]
)


def require(ok, message):
    if not ok:
        raise ValueError(message)


def clenshaw_with_derivative(coefficients, x):
    """Simultaneous backward evaluation of sum c[k]T_k(x) and its derivative."""
    b1 = b2 = d1 = d2 = 0.0
    for coefficient in reversed(coefficients[1:]):
        b = coefficient + 2 * x * b1 - b2
        d = 2 * b1 + 2 * x * d1 - d2
        b2, b1 = b1, b
        d2, d1 = d1, d
    return coefficients[0] + x * b1 - b2, b1 + x * d1 - d2


class Ephemeris:
    def __init__(self, path):
        self.path = Path(path)
        self._decode(self.path.read_bytes())

    def _decode(self, data):
        """Internal byte decoder also exercised with manufactured corrupt buffers."""
        require(
            len(data) >= 3 * RECORD_BYTES and len(data) % RECORD_BYTES == 0,
            "Incomplete coefficient record",
        )
        # Exact old x86 C-structure layout of this file, including four-byte
        # padding before timeData and AU. Not a generic SPK/DAF file reader.
        self.labels = [data[i : i + 84].split(b"\0")[0].decode("ascii") for i in (0, 84, 168)]
        self.header_start, self.header_end, self.span = struct.unpack_from("<3d", data, 2656)
        self.nconst = struct.unpack_from("<i", data, 2680)[0]
        self.file_AU_km, self.header_emrat = struct.unpack_from("<2d", data, 2688)
        self.pointers = [struct.unpack_from("<3i", data, 2704 + 12 * i) for i in range(12)]
        self.denum = struct.unpack_from("<i", data, 2848)[0]
        self.libration_pointer = struct.unpack_from("<3i", data, 2852)
        require(
            self.denum == 405 and 0 < self.nconst <= 400, "Unsupported header identity/constants"
        )
        require(
            all(
                math.isfinite(q)
                for q in (
                    self.header_start,
                    self.header_end,
                    self.span,
                    self.file_AU_km,
                    self.header_emrat,
                )
            )
            and self.header_start < self.header_end
            and self.span == 32
            and self.file_AU_km > 0
            and self.header_emrat > 0,
            "Invalid header values",
        )
        names = [
            data[252 + 6 * i : 258 + 6 * i].decode("ascii").rstrip(" \0")
            for i in range(self.nconst)
        ]
        values = struct.unpack_from("<400d", data, RECORD_BYTES)
        require(len(set(names)) == self.nconst and all(names), "Repeated/empty constant name")
        require(all(math.isfinite(v) for v in values[: self.nconst]), "Nonfinite constant")
        self.constants = dict(zip(names, values))
        require(
            self.constants["AU"] == self.file_AU_km
            and self.constants["EMRAT"] == self.header_emrat
            and self.constants["DENUM"] == 405
            and self.file_AU_km * 1000 == AU_M,
            "Constant/header mismatch",
        )
        occupied = set()
        for body, (entry, n, granules) in enumerate(self.pointers + [self.libration_pointer]):
            components = 2 if body == 11 else 3
            require(
                entry >= 3 and 1 <= n <= 50 and 1 <= granules <= 32, "Invalid coefficient pointer"
            )
            slots = set(range(entry - 1, entry - 1 + components * n * granules))
            require(
                max(slots) < 1018 and not (slots & occupied), "Coefficient pointer overlap/range"
            )
            occupied |= slots
        require(occupied == set(range(2, 1018)), "Missing coefficient coverage")
        self.records = np.frombuffer(data, dtype="<f8", offset=2 * RECORD_BYTES).reshape(-1, 1018)
        require(np.isfinite(self.records).all(), "Nonfinite coefficient record")
        starts = self.records[:, 0]
        ends = self.records[:, 1]
        require(
            np.all(ends - starts == self.span) and np.all(starts[1:] == ends[:-1]),
            "Noncontiguous/incorrect record dates",
        )
        require(
            starts[0] == self.header_start and ends[-1] <= self.header_end,
            "Record/header coverage inconsistent",
        )
        self.actual_start = float(starts[0])
        self.actual_end = float(ends[-1])
        self.metadata = dict(
            file_bytes=len(data),
            record_bytes=RECORD_BYTES,
            coefficient_records=len(self.records),
            header_JD=[self.header_start, self.header_end],
            actual_JD=[self.actual_start, self.actual_end],
            span_days=self.span,
            header_extends_past_file_days=self.header_end - self.actual_end,
            missing_advertised_records=int((self.header_end - self.actual_end) / self.span),
            AU_km=self.file_AU_km,
            EMRAT=self.header_emrat,
            constant_count=self.nconst,
            coefficient_pointers=self.pointers,
            libration_pointer=self.libration_pointer,
            full_record_finiteness_and_contiguity_verified=True,
        )

    def record_index(self, jd):
        require(
            isinstance(jd, (int, float, np.floating)) and math.isfinite(jd),
            "Nonfinite/invalid epoch",
        )
        require(
            self.actual_start <= jd <= self.actual_end
            and self.header_start <= jd <= self.header_end,
            "Epoch outside actual coefficient coverage",
        )
        return min(int(math.floor((jd - self.actual_start) / self.span)), len(self.records) - 1)

    def raw_on_record(self, body, jd, record_index, granule=None, numpy_reference=False):
        require(type(body) is int and 0 <= body <= 10, "Unsupported raw body")
        require(
            type(record_index) is int and 0 <= record_index < len(self.records),
            "Invalid record index",
        )
        record = self.records[record_index]
        start, end = record[:2]
        require(start <= jd <= end, "Epoch outside specified record")
        entry, n, groups = self.pointers[body]
        width = (end - start) / groups
        if granule is None:
            granule = min(int(math.floor((jd - start) / width)), groups - 1)
        require(type(granule) is int and 0 <= granule < groups, "Invalid granule")
        segment = start + granule * width
        require(segment <= jd <= segment + width, "Epoch outside specified granule")
        tau = 2.0 * (jd - segment) / width - 1.0
        offset = entry - 1 + 3 * n * granule
        coeff = record[offset : offset + 3 * n].reshape(3, n)
        p = []
        v = []
        for axis in coeff:
            if numpy_reference:
                pos = float(np.polynomial.chebyshev.chebval(tau, axis))
                der = float(
                    np.polynomial.chebyshev.chebval(tau, np.polynomial.chebyshev.chebder(axis))
                )
            else:
                pos, der = clenshaw_with_derivative(axis, tau)
            p.append(pos)
            v.append(der * 2.0 / width / DAY_S)
        result = np.array(p + v)
        require(np.isfinite(result).all(), "Nonfinite evaluated state")
        return result

    def raw(self, body, jd, numpy_reference=False):
        return self.raw_on_record(
            body, float(jd), self.record_index(jd), numpy_reference=numpy_reference
        )

    def heliocentric(self, body, mjd, numpy_reference=False):
        require(body in ("earth", "mercury"), "Unsupported physical body")
        require(isinstance(mjd, (int, float, np.floating)) and math.isfinite(mjd), "Invalid MJD")
        jd = float(mjd) + JD0
        target = self.raw(2 if body == "earth" else 0, jd, numpy_reference)
        if body == "earth":
            target = target - self.raw(9, jd, numpy_reference) / (1 + self.constants["EMRAT"])
        target = target - self.raw(10, jd, numpy_reference)
        return np.r_[ROT @ target[:3] / (AU_M / 1000), ROT @ target[3:] * (DAY_S / (AU_M / 1000))]

    def diagnostic_position_derivative(self, body, mjd, order):
        """NumPy higher position derivative, used only to check the FD envelope.

        Returned m/s^order. Sampling is not proof of an interval-wide bound.
        """
        require(body in ("earth", "mercury"), "Unsupported physical body")
        require(order in (2, 3), "Diagnostic acceleration/jerk only")
        jd = float(mjd) + JD0
        ri = self.record_index(jd)
        record = self.records[ri]

        def raw_derivative(index):
            entry, n, groups = self.pointers[index]
            width = self.span / groups
            group = min(int((jd - record[0]) / width), groups - 1)
            tau = 2 * (jd - (record[0] + group * width)) / width - 1
            offset = entry - 1 + 3 * n * group
            coeff = record[offset : offset + 3 * n].reshape(3, n)
            return (
                np.array(
                    [
                        np.polynomial.chebyshev.chebval(
                            tau, np.polynomial.chebyshev.chebder(c, m=order)
                        )
                        for c in coeff
                    ]
                )
                * (2 / width / DAY_S) ** order
            )

        derivative = raw_derivative(2 if body == "earth" else 0) - raw_derivative(10)
        if body == "earth":
            derivative -= raw_derivative(9) / (1 + self.header_emrat)
        return ROT @ derivative * 1000


DEPARTURE_MJD = 52654.5


DURATION_DOMAIN = (0.5, 600.0)


RAW_BODIES = (0, 10)


def _scalar(value, label):
    require(
        isinstance(value, (int, float, np.integer, np.floating))
        and not isinstance(value, (bool, np.bool_)),
        label + " must be a real scalar",
    )
    result = float(value)
    require(math.isfinite(result), label + " must be finite")
    return result


def _departure(value):
    result = _scalar(value, "Departure MJD")
    require(result == DEPARTURE_MJD, "Departure is frozen for this Mercury development model")
    return result


def _duration_range(value):
    require(
        isinstance(value, (tuple, list)) and len(value) == 2,
        "Duration range must contain exactly two endpoints",
    )
    left, right = (_scalar(v, "Duration endpoint") for v in value)
    require(
        DURATION_DOMAIN[0] <= left < right <= DURATION_DOMAIN[1],
        "Duration range lies outside the declared development domain",
    )
    return left, right


def _load_ephemeris(data_path):
    eph = Ephemeris(Path(data_path))
    require(ca.__version__ == "3.7.2", "Model is pinned to CasADi 3.7.2")
    return eph


def _coefficient_pieces(eph, left, right):
    """Choose coefficients by a strict interior midpoint, never by a join tie."""
    mid_jd = left + (right - left) / 2.0 + JD0
    record_index = eph.record_index(mid_jd)
    record = eph.records[record_index]
    result = []
    for body in RAW_BODIES:
        entry, count, groups = eph.pointers[body]
        width = float((record[1] - record[0]) / groups)
        group = int(math.floor((mid_jd - record[0]) / width))
        require(0 <= group < groups, "Coefficient midpoint group out of range")
        start = float(record[0] - JD0 + group * width)
        require(
            start <= left < right <= start + width,
            "Selected common cell crosses a raw coefficient interval",
        )
        offset = entry - 1 + 3 * count * group
        coefficients = record[offset : offset + 3 * count].reshape(3, count).copy()
        coefficients.flags.writeable = False
        result.append(
            dict(
                body=body,
                group=group,
                width=width,
                start_MJD=start,
                end_MJD=start + width,
                coefficients=coefficients,
                coefficients_per_axis=count,
            )
        )
    require(
        max(p["start_MJD"] for p in result) == left and min(p["end_MJD"] for p in result) == right,
        "Bounds are not one exact common Mercury/Sun coefficient cell",
    )
    return record_index, result


def _coefficient_descriptors(eph, departure, duration_range):
    lo, hi = duration_range
    # Work with exact half-integer MJD/development endpoints only for finding
    # coefficients. Evaluation below uses small duration-relative offsets.
    epoch_lo, epoch_hi = departure + lo, departure + hi
    require(
        eph.actual_start - JD0 <= epoch_lo < epoch_hi <= eph.actual_end - JD0,
        "Requested span exceeds actual file coverage",
    )
    edges = set()
    for record in eph.records:
        start, end = (float(q - JD0) for q in record[:2])
        if end < epoch_lo or start > epoch_hi:
            continue
        for body in RAW_BODIES:
            groups = eph.pointers[body][2]
            width = (end - start) / groups
            edges.update(start + k * width for k in range(groups + 1))
    ordered = sorted(edges)
    cells = []
    for left, right in zip(ordered, ordered[1:]):
        full_lo, full_hi = left - departure, right - departure
        bounded_lo, bounded_hi = max(lo, full_lo), min(hi, full_hi)
        if bounded_lo >= bounded_hi:
            continue  # Do not invent zero-width cells at the requested endpoints.
        record_index, pieces = _coefficient_pieces(eph, left, right)
        cells.append(
            dict(
                schema=1,
                target="mercury",
                departure_MJD=departure,
                requested_duration_range_days=[lo, hi],
                cell_MJD=[left, right],
                coefficient_duration_bounds_days=[full_lo, full_hi],
                duration_bounds_days=[bounded_lo, bounded_hi],
                record_index=record_index,
                raw_body_indices=list(RAW_BODIES),
                coefficient_intervals=[
                    {
                        k: p[k]
                        for k in ("body", "group", "start_MJD", "end_MJD", "coefficients_per_axis")
                    }
                    for p in pieces
                ],
            )
        )
    require(
        bool(cells)
        and cells[0]["duration_bounds_days"][0] == lo
        and cells[-1]["duration_bounds_days"][1] == hi,
        "Incomplete requested cell coverage",
    )
    require(
        all(
            a["duration_bounds_days"][1] == b["duration_bounds_days"][0]
            for a, b in zip(cells, cells[1:])
        ),
        "Noncontiguous common-cell coverage",
    )
    return cells


def enumerate_cells(data_path, departure=DEPARTURE_MJD, duration_range=DURATION_DOMAIN):
    """Return ordered JSON-safe exact-cell descriptors covering a duration range.

    Full coefficient bounds and the clipped admissible duration interval are
    separate fields. The latter, not the full cell, is the permitted solver
    domain. Adjacent closed domains overlap at their shared endpoint only.
    """
    departure = _departure(departure)
    duration_range = _duration_range(duration_range)
    eph = _load_ephemeris(data_path)
    return _coefficient_descriptors(eph, departure, duration_range)


def select_cell(cells, duration, *, side):
    """Explicit numeric side selection; this operation is not differentiated.

    At an interior join, ``left`` chooses the preceding polynomial and ``right``
    the following one. Away from a join (including the requested outer bounds),
    the only available containing cell is returned. A list that was modified
    after enumeration is still rejected by ``build_cell`` before evaluation.
    """
    require(side in ("left", "right"), "Join side must explicitly be left or right")
    require(
        isinstance(cells, (tuple, list)) and len(cells) > 0,
        "Expected a nonempty list of cell descriptors",
    )
    duration = _scalar(duration, "Duration")
    bounds = []
    for cell in cells:
        require(isinstance(cell, dict), "Malformed cell descriptor")
        bound = cell.get("duration_bounds_days")
        require(isinstance(bound, list) and len(bound) == 2, "Malformed cell bounds")
        left, right = (_scalar(v, "Cell bound") for v in bound)
        require(left < right, "Empty or reversed cell domain")
        bounds.append((left, right))
    require(
        all(a[1] == b[0] for a, b in zip(bounds, bounds[1:])),
        "Cell list must be ordered and contiguous",
    )
    indices = [i for i, (left, right) in enumerate(bounds) if left <= duration <= right]
    require(len(indices) in (1, 2), "Duration outside selected cell coverage")
    return cells[indices[0] if side == "left" else indices[-1]]


def build_cell(data_path, cell, departure=DEPARTURE_MJD):
    """Build one immutable-in-use coefficient graph and bounded NumPy wrappers.

    ``state``, ``derivative`` and ``second_derivative`` are local CasADi Functions
    of T_days. Their raw mathematical expressions do NOT enforce bounds; NLP
    callers must apply returned bounds. All ``evaluate*`` and ``independent*``
    numeric wrappers enforce them. The independent jets use differentiated
    NumPy Chebyshev coefficients rather than differentiating the CasADi graph.
    """
    departure = _departure(departure)
    require(isinstance(cell, dict), "Cell must be an enumerated descriptor")
    requested_range = _duration_range(cell.get("requested_duration_range_days"))
    eph = _load_ephemeris(data_path)
    candidates = _coefficient_descriptors(eph, departure, requested_range)
    # Dict equality by itself accepts True == 1. Require JSON field types as
    # well as exact values by comparing canonical serialisations.
    import json

    try:
        serial = json.dumps(cell, sort_keys=True, allow_nan=False)
    except (TypeError, ValueError) as error:
        raise ValueError("Cell descriptor must be finite JSON data") from error
    matching = [d for d in candidates if json.dumps(d, sort_keys=True, allow_nan=False) == serial]
    require(len(matching) == 1, "Cell descriptor differs from the pinned exact coefficient map")
    descriptor = matching[0]
    left, right = descriptor["cell_MJD"]
    bounds = tuple(descriptor["duration_bounds_days"])
    full_bounds = tuple(descriptor["coefficient_duration_bounds_days"])
    record_index, pieces = _coefficient_pieces(eph, left, right)
    duration = ca.SX.sym("duration_days")
    local_offset = duration - full_bounds[0]

    def symbolic_position(piece):
        tau = 2.0 * ((left - piece["start_MJD"]) + local_offset) / piece["width"] - 1.0
        axes = []
        for coefficients in piece["coefficients"]:
            t0, t1 = ca.SX(1), tau
            value = float(coefficients[0]) + float(coefficients[1]) * t1
            for k in range(2, len(coefficients)):
                t2 = 2.0 * tau * t1 - t0
                value += float(coefficients[k]) * t2
                t0, t1 = t1, t2
            axes.append(value)
        return ca.vertcat(*axes)

    position = ca.mtimes(
        ca.DM(ROT), symbolic_position(pieces[0]) - symbolic_position(pieces[1])
    ) / (AU_M / 1000.0)
    velocity = ca.jacobian(position, duration)
    state_expr = ca.vertcat(position, velocity)
    derivative_expr = ca.jacobian(state_expr, duration)
    second_expr = ca.jacobian(derivative_expr, duration)
    name = "mercury_cell_" + str(record_index) + "_" + str(pieces[0]["group"])
    state = ca.Function(name + "_state", [duration], [state_expr])
    derivative = ca.Function(name + "_derivative", [duration], [derivative_expr])
    second_derivative = ca.Function(name + "_second_derivative", [duration], [second_expr])

    def checked_duration(value):
        value = _scalar(value, "Duration")
        require(bounds[0] <= value <= bounds[1], "Duration outside this closed admissible cell")
        return value

    def checked_array(value):
        result = np.asarray(value, dtype=float).reshape(6)
        require(np.isfinite(result).all(), "Nonfinite target state or derivative")
        return result

    def evaluate(value):
        return checked_array(state(checked_duration(value)))

    def evaluate_derivative(value):
        return checked_array(derivative(checked_duration(value)))

    def evaluate_second_derivative(value):
        return checked_array(second_derivative(checked_duration(value)))

    def independent_orders(value):
        """Return r, v, a, jerk in AU/day**order at the same smooth local epoch."""
        value = checked_duration(value)
        sums = [np.zeros(3) for _ in range(4)]
        for sign, piece in zip((1.0, -1.0), pieces):
            tau = (
                2.0 * ((left - piece["start_MJD"]) + (value - full_bounds[0])) / piece["width"]
                - 1.0
            )
            for order in range(4):
                raw = (
                    np.array(
                        [
                            np.polynomial.chebyshev.chebval(
                                tau, np.polynomial.chebyshev.chebder(axis, m=order)
                            )
                            for axis in piece["coefficients"]
                        ]
                    )
                    * (2.0 / piece["width"]) ** order
                )
                sums[order] += sign * raw
        result = [ROT @ value / (AU_M / 1000.0) for value in sums]
        require(all(np.isfinite(value).all() for value in result), "Nonfinite independent jet")
        return result

    def independent_state(value):
        orders = independent_orders(value)
        return checked_array(np.r_[orders[0], orders[1]])

    def independent_derivative(value):
        orders = independent_orders(value)
        return checked_array(np.r_[orders[1], orders[2]])

    def independent_second_derivative(value):
        orders = independent_orders(value)
        return checked_array(np.r_[orders[2], orders[3]])

    def native_state(value):
        """Existing decoder with native MJD->JD rounding/right-side dispatch.

        This is NOT necessarily the selected cell's left limit at a join and is
        NOT necessarily the same smooth epoch. Align epoch and side explicitly
        when treating its difference as arithmetic/model error.
        """
        return checked_array(eph.heliocentric("mercury", departure + checked_duration(value)))

    # Triangle inequality and |Tk(tau)| <= 1 provide coefficient-based upper
    # bounds over the FULL selected raw cells, not sample-derived assertions.
    derivative_bounds_SI = {}
    for order in (1, 2, 3, 4, 5):
        bound = 0.0
        for piece in pieces:
            axis_bounds = np.array(
                [
                    np.sum(np.abs(np.polynomial.chebyshev.chebder(axis, m=order)))
                    for axis in piece["coefficients"]
                ]
            )
            bound += float(np.linalg.norm(axis_bounds)) * (2.0 / piece["width"]) ** order
        derivative_bounds_SI[str(order)] = bound * 1000.0 / DAY_S**order

    provenance = dict(
        descriptor,
        ephemeris="DE405",
        CasADi_version=ca.__version__,
        state_units=["AU"] * 3 + ["AU/day"] * 3,
        derivative_units=["AU/day"] * 3 + ["AU/day^2"] * 3,
        second_derivative_units=["AU/day^2"] * 3 + ["AU/day^3"] * 3,
        position_derivative_bounds_SI=derivative_bounds_SI,
        bound_derivation=(
            "Sum of per-body Euclidean norms of per-axis absolute differentiated Chebyshev "
            "coefficient sums; |Tk|<=1 and orthogonal frame rotation"
        ),
        duration_offset_rule=(
            "(T - full_cell_lower_T) + (cell_left_MJD - raw_piece_start_MJD); "
            "no large JD in symbolic evaluation"
        ),
        numeric_wrappers_enforce_bounds=True,
        symbolic_bounds_enforced_by_caller=True,
        cell_endpoints="closed, one-sided polynomial values at coefficient joins",
        global_differentiability_claim=False,
        global_dispatch_differentiated=False,
        native_execution=False,
        optimiser_execution=False,
    )
    return dict(
        state=state,
        derivative=derivative,
        second_derivative=second_derivative,
        evaluate=evaluate,
        evaluate_derivative=evaluate_derivative,
        evaluate_second_derivative=evaluate_second_derivative,
        independent=independent_state,
        independent_state=independent_state,
        independent_derivative=independent_derivative,
        independent_second_derivative=independent_second_derivative,
        independent_orders=independent_orders,
        native_state=native_state,
        cell=tuple(descriptor["cell_MJD"]),
        duration_bounds_days=bounds,
        coefficient_duration_bounds_days=full_bounds,
        provenance=provenance,
    )


SCHEMA = "fresh600_declared_time_domain_v2"


CELLS = ((53248.0, 53256.0), (53240.0, 53248.0), (53232.0, 53240.0), (53224.0, 53232.0))


DURATIONS = ((593.5, 600.0), (585.5, 593.5), (577.5, 585.5), (569.5, 577.5))


GRID_ANCHOR_MJD = 53248.0


GRID_WIDTH_DAYS = 8.0


def _pair(value, name):
    if not isinstance(value, (tuple, list)) or len(value) != 2:
        raise ValueError(name + " must contain two declared numeric bounds")
    try:
        valid = all(type(x) in (int, float) and math.isfinite(x) for x in value)
    except OverflowError:
        valid = False
    if not valid:
        raise ValueError(name + " must contain finite numbers, not booleans")
    return tuple(float(x) for x in value)


def validate_domain(continuation_domain=None):
    if continuation_domain is None:
        return None
    fields = {"minimum_duration_days", "maximum_duration_days"}
    if type(continuation_domain) is not dict or set(continuation_domain) != fields:
        raise ValueError("Continuation domain requires exactly minimum/maximum duration")
    lo, hi = _pair(
        [
            continuation_domain["minimum_duration_days"],
            continuation_domain["maximum_duration_days"],
        ],
        "continuation domain",
    )
    model_lo, model_hi = DURATION_DOMAIN
    if not model_lo <= lo < hi <= model_hi:
        raise ValueError("Continuation domain lies outside the existing0.5..600-day model")
    return {"minimum_duration_days": lo, "maximum_duration_days": hi}


def _grid(continuation_domain):
    # Check the declaration once, then use the same range throughout this call.
    contract = validate_domain(continuation_domain)
    lo, hi = (
        (569.5, 600.0)
        if contract is None
        else (contract["minimum_duration_days"], contract["maximum_duration_days"])
    )
    epoch_lo, epoch_hi = (DEPARTURE_MJD + value for value in (lo, hi))
    index = math.ceil((epoch_hi - GRID_ANCHOR_MJD) / GRID_WIDTH_DAYS) - 1
    cells = []
    while True:
        left = GRID_ANCHOR_MJD + index * GRID_WIDTH_DAYS
        right = left + GRID_WIDTH_DAYS
        if right <= epoch_lo:
            break
        if left < epoch_hi:
            cells.append((left, right))
        index -= 1
    if not cells:
        raise ValueError("Empty continuation cell domain")
    return contract, lo, hi, tuple(cells)


def _domain_record(contract, lo: float, hi: float, cell, index: int):
    raw = tuple(value - DEPARTURE_MJD for value in cell)
    return {
        "schema": SCHEMA,
        "index": index,
        "cell_MJD": list(cell),
        "duration_bounds_days": [max(lo, raw[0]), min(hi, raw[1])],
        "coefficient_duration_bounds_days": list(raw),
        "continuation_domain": copy.deepcopy(contract),
    }


def campaign_cells(continuation_domain=None):
    """Return high-to-low raw coefficient cells intersecting the declared range."""
    return _grid(continuation_domain)[3]


def campaign_domains(continuation_domain=None):
    """Describe all cells without rebuilding the same grid for every cell."""
    contract, lo, hi, cells = _grid(continuation_domain)
    return [_domain_record(contract, lo, hi, cell, index) for index, cell in enumerate(cells)]


def domain(cell_MJD=None, continuation_domain=None):
    contract, lo, hi, cells = _grid(continuation_domain)
    cell = cells[0] if cell_MJD is None else _pair(cell_MJD, "arrival cell")
    if cell not in cells:
        raise ValueError("Arrival cell is outside the declared Mercury domains")
    return _domain_record(contract, lo, hi, cell, cells.index(cell))


def bounds(candidate):
    """Earlier/explicit-domain candidates must carry their actual clipped bounds."""
    if type(candidate) is not dict or "arrival_cell_MJD" not in candidate:
        raise ValueError("Candidate must identify its exact arrival cell")
    selected = domain(candidate["arrival_cell_MJD"], candidate.get("continuation_domain"))
    expected = tuple(selected["duration_bounds_days"])
    if "duration_bounds_days" not in candidate:
        if selected["index"] != 0 or selected["continuation_domain"] is not None:
            raise ValueError(
                "Earlier/declared-domain candidate requires explicit duration_bounds_days"
            )
    elif _pair(candidate["duration_bounds_days"], "duration bounds") != expected:
        raise ValueError("Candidate duration bounds differ from its exact cell domain")
    return expected


def descriptors(ephemeris_path, continuation_domain=None):
    """Match the full declared grid and clipped bounds to actual coefficient data."""
    return _campaign_descriptors(ephemeris_path, _grid(continuation_domain))


def _campaign_descriptors(ephemeris_path, grid):
    contract, lo, hi, cells = grid
    actual = enumerate_cells(ephemeris_path, duration_range=(lo, hi))
    expected = [_domain_record(contract, lo, hi, cell, index) for index, cell in enumerate(cells)]
    if len(actual) != len(expected):
        raise ValueError("DE405 common-cell count differs from the declared contract")
    ordered = []
    for selected in expected:
        matches = [
            item
            for item in actual
            if item["cell_MJD"] == selected["cell_MJD"]
            and item["duration_bounds_days"] == selected["duration_bounds_days"]
            and item["departure_MJD"] == DEPARTURE_MJD
            and item["target"] == "mercury"
            and item["raw_body_indices"] == [0, 10]
        ]
        if len(matches) != 1:
            raise ValueError("DE405 metadata differs from a declared Mercury/Sun cell")
        ordered.append(copy.deepcopy(matches[0]))
    return ordered


def descriptor(ephemeris_path, cell_MJD=None, continuation_domain=None):
    grid = _grid(continuation_domain)
    cells = grid[3]
    cell = cells[0] if cell_MJD is None else _pair(cell_MJD, "arrival cell")
    if cell not in cells:
        raise ValueError("Arrival cell is outside the declared Mercury domains")
    return _campaign_descriptors(ephemeris_path, grid)[cells.index(cell)]


def neighbours(cell_MJD=None, continuation_domain=None):
    cells = _grid(continuation_domain)[3]
    selected = cells[0] if cell_MJD is None else _pair(cell_MJD, "arrival cell")
    if selected not in cells:
        raise ValueError("Arrival cell is outside the declared Mercury domains")
    return [
        list(cell)
        for cell in cells
        if cell != selected and (cell[1] == selected[0] or selected[1] == cell[0])
    ]
