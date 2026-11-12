"""Keep a checked point; move to the earlier target cell only when needed."""

from __future__ import annotations
import copy
import hashlib
import json
import math
from pathlib import Path
import time
import numpy as np
from . import records
from .records import require


EXPECTED_REJECTIONS = {
    "Maximum_WallTime_Exceeded",
    "Maximum_CpuTime_Exceeded",
    "Maximum_Iterations_Exceeded",
    "Infeasible_Problem_Detected",
    "Restoration_Failed",
    "Search_Direction_Becomes_Too_Small",
}
TIME_TIE_DAYS = 1e-6
BUDGET_STATUSES = {
    "Maximum_WallTime_Exceeded",
    "Maximum_CpuTime_Exceeded",
    "Maximum_Iterations_Exceeded",
}
ADMISSION_POLICY = "converged_or_budget_feasible_v2"
MINIMUM_ALLOCATION_SECONDS = 180.0
NON_SOLVER_ALLOWANCE_SECONDS = 120.0


HANDOVER_SCHEMA = "explicit_earlier_cell_handover_v1"
MAX_RETIME_DAYS = 1e-6


def solver_admissible(result):
    """A resource stop can return a feasible point, not an optimality certificate."""
    if result.get("solver_converged") is True and result.get("solver_status") == "Solve_Succeeded":
        return True
    return bool(
        result.get("admission_policy") == ADMISSION_POLICY
        and result.get("phase") == "quality"
        and result.get("solver_status") in BUDGET_STATUSES
        and result.get("complete_convergence_metadata") is True
        and result.get("solver_converged") is False
        and result.get("solver_reported_success") is False
        and result.get("finite_returned_vector") is True
        and result.get("objective_not_increased") is True
        and not result.get("errors")
    )


def admitted(bundle):
    if not isinstance(bundle, dict):
        return False
    result, host, candidate = (
        bundle.get("result", {}),
        bundle.get("host_validation") or {},
        bundle.get("candidate"),
    )
    violation = result.get("constraint_violation")
    return bool(
        result.get("eligible_for_native_review") is True
        and solver_admissible(result)
        and result.get("original_bounds_met") is True
        and type(violation) in (int, float)
        and math.isfinite(violation)
        and 0 <= violation <= 1e-12
        and host.get("accepted_for_native_review") is True
        and isinstance(candidate, dict)
    )


def _quality(bundle):
    candidate = bundle["candidate"]
    goal = candidate["mission_profile"]
    encounter = bundle["host_validation"]["encounter"]
    values = (
        encounter["position_miss_km"] * 1000.0 / goal["position_target_m"],
        encounter["velocity_miss_m_s"] / goal["velocity_target_m_s"],
    )
    require(
        all(
            type(value) in (int, float) and math.isfinite(value) and value >= 0 for value in values
        ),
        "Invalid candidate quality metrics",
    )
    return max(values)


def choose(restoration, refinement):
    """Rank host-admitted candidates; final native replay is still mandatory."""
    require(admitted(restoration), "No admitted same-run incumbent")
    require(
        isinstance(refinement, dict)
        and type(refinement.get("result", {}).get("solver_calls")) is int
        and refinement["result"]["solver_calls"] == 1,
        "Refinement is not a recorded single solver call",
    )
    incumbent = restoration["candidate"]
    trial = refinement.get("candidate")
    if trial is not None:
        require(
            all(
                trial.get(key) == incumbent.get(key)
                for key in (
                    "run_id",
                    "origin_controller_sha256",
                    "mission_profile",
                    "local_protocol_sha256",
                )
            ),
            "Refinement and incumbent belong to different runs/contracts",
        )
    selected, reason = "restore", "expected_refinement_rejection"
    refinement_admitted = admitted(refinement)
    if refinement_admitted:
        difference = trial["duration_days"] - incumbent["duration_days"]
        if difference < -TIME_TIE_DAYS:
            selected, reason = "quality", "shorter_at_declared_quality"
        elif abs(difference) <= TIME_TIE_DAYS and _quality(refinement) < _quality(restoration):
            selected, reason = "quality", "equal_time_lower_normalized_error"
        else:
            reason = "refinement_did_not_improve_declared_order"
    else:
        result = refinement.get("result", {})
        require(
            result.get("complete_convergence_metadata") is True
            and result.get("solver_status") in EXPECTED_REJECTIONS
            and not result.get("errors"),
            "Unexpected solver/validation failure requires investigation",
        )
    chosen = refinement if selected == "quality" else restoration
    return copy.deepcopy(chosen["candidate"]), {
        "schema": "same_run_incumbent_selection_v2",
        "selected_phase": selected,
        "reason": reason,
        "refinement_status": refinement["result"]["solver_status"],
        "refinement_admitted": refinement_admitted,
        "incumbent_admitted": True,
        "historical_fallback": False,
        "native_validation_pending": True,
        "selected_solver_converged": chosen["result"]["solver_converged"],
        "selected_admission": (
            "converged" if chosen["result"]["solver_converged"] else "budget_feasible_not_converged"
        ),
        "local_or_global_optimality_established": False,
        "time_tie_days": TIME_TIE_DAYS,
    }


class BudgetExpired(RuntimeError):
    """Only a recorded resource stop, never a disguised scientific exception."""

    def __init__(self, detail=None):
        super().__init__("Declared optimisation wall budget exhausted")
        self.detail = detail or {}


class Budget:
    # The clock is a brake. We do not have to use all 50 minutes.
    def __init__(self, seconds=3000.0, clock=time.monotonic):
        require(
            type(seconds) in (int, float) and math.isfinite(seconds) and seconds > 0,
            "Invalid optimisation wall budget",
        )
        self.clock, self.limit, self.start = clock, float(seconds), clock()
        self.stopped = None

    def elapsed(self):
        return max(0.0, (self.clock() if self.stopped is None else self.stopped) - self.start)

    def stop(self):
        if self.stopped is None:
            self.stopped = self.clock()

    def remaining(self):
        return max(0.0, self.limit - self.elapsed())

    def allocation(self):
        remaining = self.remaining()
        # Includes graph building, validation and saving inside the outer budget.
        if remaining < MINIMUM_ALLOCATION_SECONDS:
            raise BudgetExpired()
        return {
            "solver_seconds": min(900.0, remaining - NON_SOLVER_ALLOWANCE_SECONDS),
            "worker_seconds": min(1200.0, remaining),
            "remaining_seconds": remaining,
        }


def key(index, attempt):
    require(
        type(index) is int and index >= 0 and type(attempt) is int and 0 <= attempt < 3,
        "Invalid cell/attempt index",
    )
    return f"quality_{index + 1:03d}_{attempt + 1:02d}"


def run_with_budget(restoration, cells, solve, handover, checkpoint, budget):
    require(admitted(restoration), "No guarded same-run restoration")
    incumbent, selected_key = restoration, "restore"
    current = copy.deepcopy(restoration["candidate"])
    ledger, last_status, stop = [], "NOT_ATTEMPTED", "declared_domain_exhausted"
    checkpoint(incumbent, selected_key, ledger, "restoration_admitted")
    finished = False
    for index, cell in enumerate(cells):
        require(current["arrival_cell_MJD"] == cell, "Continuation input cell mismatch")
        for attempt in range(3):
            try:
                allocation = budget.allocation()
                trial = solve(index, attempt, current, allocation)
            except BudgetExpired as error:
                if error.detail:
                    ledger.append(dict(error.detail, cell_index=index, attempt_index=attempt))
                stop, finished = (
                    error.detail.get("stop_reason", "optimisation_wall_budget_exhausted"),
                    True,
                )
                break
            require(
                trial.get("candidate") is None or trial["candidate"]["arrival_cell_MJD"] == cell,
                "Solver returned another cell",
            )
            last_status = trial["result"]["solver_status"]
            previous_time = incumbent["candidate"]["duration_days"]
            chosen, decision = choose(incumbent, trial)
            improved = decision["selected_phase"] == "quality"
            if improved:
                incumbent, selected_key = trial, key(index, attempt)
            entry = {
                "stage_key": key(index, attempt),
                "cell_index": index,
                "attempt_index": attempt,
                "solver_calls": trial["result"]["solver_calls"],
                "solver_status": last_status,
                "trial_admitted": decision["refinement_admitted"],
                "selected_duration_days": chosen["duration_days"],
                "selected_stage_key": selected_key,
                "allocation": allocation,
                "selection": decision,
            }
            ledger.append(entry)
            checkpoint(incumbent, selected_key, ledger, "attempt_recorded")
            if not entry["trial_admitted"] or not improved:
                stop, finished = "no_better_admitted_candidate", True
                break
            floor = chosen["duration_bounds_days"][0]
            if 0.0 <= chosen["duration_days"] - floor <= 1e-6:
                if index == len(cells) - 1:
                    finished = True
                    break
                remaining = budget.remaining()
                if remaining < MINIMUM_ALLOCATION_SECONDS:
                    # A handover cannot make an already unaffordable next solve
                    # affordable. Retain this admitted point without doing work
                    # whose result cannot be used under the unchanged budget.
                    stop, finished = "insufficient_local_budget_before_handover", True
                    entry["handover_skipped"] = {
                        "reason": stop,
                        "remaining_seconds": remaining,
                        "minimum_next_solve_allocation_seconds": MINIMUM_ALLOCATION_SECONDS,
                    }
                    break
                moved = handover(index + 1, chosen, cells[index + 1], selected_key)
                entry["handover"] = moved["record"]
                checkpoint(incumbent, selected_key, ledger, "handover_recorded")
                if (moved.get("host_validation") or {}).get(
                    "accepted_for_native_review"
                ) is not True:
                    stop, finished = moved["record"].get("reason", "handover_rejected"), True
                    break
                require(
                    all(
                        moved["candidate"][field] == chosen[field]
                        for field in (
                            "run_id",
                            "origin_controller_sha256",
                            "mission_profile",
                            "local_protocol_sha256",
                        )
                    ),
                    "Handover changed campaign identity",
                )
                current = moved["candidate"]
                break
            # An interior converged point ends this start; more time alone is not a new trial.
            if trial["result"]["solver_converged"]:
                stop, finished = "converged_interior_local_solution", True
                break
            # Retry only a progressing, physically checked resource-limited point.
            require(last_status in BUDGET_STATUSES, "Unclassified continuation state")
            if chosen["duration_days"] >= previous_time - TIME_TIE_DAYS:
                stop, finished = "no_time_progress_for_budget_continuation", True
                break
            current = copy.deepcopy(chosen)
            if attempt == 2:
                stop, finished = "declared_same_cell_attempt_limit", True
        if finished:
            break
    checkpoint(incumbent, selected_key, ledger, stop)
    budget.stop()
    calls = sum(item.get("solver_calls", 0) for item in ledger)
    return copy.deepcopy(incumbent["candidate"]), {
        "schema": "long_same_run_selection_v1",
        "selected_stage_key": selected_key,
        "selected_phase": "restore" if selected_key == "restore" else "quality",
        "reason": stop,
        "refinement_status": last_status,
        "refinement_admitted": selected_key != "restore",
        "selected_solver_converged": incumbent["result"]["solver_converged"],
        "selected_admission": (
            "converged"
            if incumbent["result"]["solver_converged"]
            else "budget_feasible_not_converged"
        ),
        "quality_attempts": len(ledger),
        "quality_solver_calls": calls,
        "total_solver_calls": 1 + calls,
        "optimisation_wall_seconds": budget.elapsed(),
        "optimisation_budget_overrun_seconds": max(0.0, budget.elapsed() - budget.limit),
        "optimisation_wall_budget_seconds": budget.limit,
        "cells": ledger,
        "native_validation_pending": True,
        "historical_fallback": False,
        "local_or_global_optimality_established": False,
    }


def _digest(value):
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()
    ).hexdigest()


def _bindings(ephemeris_path):
    files = sorted(Path(__file__).parent.glob("*.py")) + [Path(ephemeris_path).resolve()]
    return {str(path.resolve()): records.sha(path) for path in files}


def _rk16_nodes(controls, duration_days):
    """Single initialization, sequential fixed-hold NLP flow; no node resets."""
    from . import model, optimiser, transfer, validation

    flow = model.build_segment(16)["flow"]
    nodes = [validation.INITIAL[:6].copy()]
    step = duration_days / optimiser.HOLDS
    for control in controls:
        following = np.asarray(
            flow(nodes[-1], control, step, transfer.LIGHTNESS), dtype=float
        ).reshape(6)
        require(np.isfinite(following).all(), "Nonfinite boundary reconstruction")
        nodes.append(following)
    values = np.asarray(nodes)
    require(
        values.shape == (optimiser.HOLDS + 1, 6)
        and np.array_equal(values[0], validation.INITIAL[:6]),
        "Boundary reconstruction initial/count mismatch",
    )
    return values


def prepare_handover(candidate, ephemeris_path, next_cell_MJD):
    """Return candidate/host_validation/record, or a recorded expected rejection.

    The parent admits the source solver result and binds this record to the
    campaign. A string on a candidate is not itself proof of earlier admission.
    Software or malformed-input errors raise; a completed failed host check is
    an expected handover rejection and must retain the incumbent.
    """
    from . import ephemeris, optimiser, validation

    started = time.monotonic()
    original = copy.deepcopy(candidate)
    profile = optimiser._identity(original)
    duration, _, controls, _ = validation.checked_candidate(
        original, "quality", mission_profile=profile
    )
    context = original.get("continuation_domain")
    source = ephemeris.domain(original["arrival_cell_MJD"], continuation_domain=context)
    destination = ephemeris.domain(next_cell_MJD, continuation_domain=context)
    cells = [tuple(cell) for cell in ephemeris.campaign_cells(context)]
    index = cells.index(tuple(source["cell_MJD"]))
    require(
        index + 1 < len(cells) and tuple(destination["cell_MJD"]) == cells[index + 1],
        "Only the immediate earlier declared cell is allowed",
    )
    boundary = source["duration_bounds_days"][0]
    require(boundary == destination["duration_bounds_days"][1], "Selected cells are not adjacent")
    offset = duration - boundary
    record = {
        "schema": HANDOVER_SCHEMA,
        "source_candidate_sha256": _digest(original),
        "source_solver_status": original.get("solver_status"),
        "source_solver_phase": original.get("solver_phase"),
        "run_id": original["run_id"],
        "origin_controller_sha256": original["origin_controller_sha256"],
        "mission_profile": profile,
        "local_protocol_sha256": original["local_protocol_sha256"],
        "continuation_domain": context,
        "control_policy": original.get("control_policy"),
        "source_cell": source,
        "destination_cell": destination,
        "source_duration_days": duration,
        "boundary_duration_days": boundary,
        "retime_days": boundary - duration,
        "retime_seconds": (boundary - duration) * 86400.0,
        "maximum_retime_days": MAX_RETIME_DAYS,
        "solver_calls": 0,
        "native_calls": 0,
        "native_mission_accepted": False,
        "source_solver_admission_owned_by_parent": True,
        "same_controls": True,
        "same_original_initial_state": True,
        "one_initialization": True,
        "shooting_resets": 0,
        "new_nodes": "sequential existing NLP RK16 segment flow; not finer-mesh shooting nodes",
    }
    if not 0.0 <= offset <= MAX_RETIME_DAYS:
        record.update(
            status="HANDOVER_REJECTED",
            reason="source_not_at_declared_cell_floor",
            accepted_for_next_cell=False,
            wall_seconds=time.monotonic() - started,
        )
        return {"candidate": None, "host_validation": None, "record": record}
    identities = _bindings(ephemeris_path)
    source_host = validation.validate_candidate(
        original, ephemeris_path, stage="quality", mission_profile=profile
    )
    if source_host.get("accepted_for_native_review") is not True:
        require(
            all(records.sha(path) == digest for path, digest in identities.items()),
            "Handover inputs changed",
        )
        record.update(
            status="HANDOVER_REJECTED",
            reason="source_host_recheck_failed",
            accepted_for_next_cell=False,
            source_host_validation=source_host,
            evidence_sha256=identities,
            wall_seconds=time.monotonic() - started,
        )
        return {"candidate": None, "host_validation": None, "record": record}
    old_arrival = optimiser.build_arrival(
        ephemeris_path, source["cell_MJD"], continuation_domain=context
    )
    new_arrival = optimiser.build_arrival(
        ephemeris_path, destination["cell_MJD"], continuation_domain=context
    )
    source_target_audit = validation.validate_target_model(
        ephemeris_path, old_arrival, continuation_domain=context
    )
    destination_target_audit = validation.validate_target_model(
        ephemeris_path, new_arrival, continuation_domain=context
    )
    old_target = validation.independent_target(
        boundary, ephemeris_path, cell_MJD=source["cell_MJD"], continuation_domain=context
    )
    new_target = validation.independent_target(
        boundary, ephemeris_path, cell_MJD=destination["cell_MJD"], continuation_domain=context
    )
    join = validation._error(old_target["state"][None, :], new_target["state"][None, :])
    join_pass = join["position_m"] <= 0.001 and join["velocity_m_s"] <= 1e-7
    if not join_pass:
        require(
            all(records.sha(path) == digest for path, digest in identities.items()),
            "Handover inputs changed",
        )
        record.update(
            status="HANDOVER_REJECTED",
            reason="target_cell_join_check_failed",
            accepted_for_next_cell=False,
            target_join=join,
            evidence_sha256=identities,
            wall_seconds=time.monotonic() - started,
        )
        return {"candidate": None, "host_validation": None, "record": record}
    # Controls stay fixed here; rebuild the nodes before using the earlier target cell.
    nodes = _rk16_nodes(controls, boundary)
    result = copy.deepcopy(original)
    result.update(
        duration_days=boundary,
        arrival_cell_MJD=destination["cell_MJD"],
        duration_bounds_days=destination["duration_bounds_days"],
        nodes_cartesian=nodes.tolist(),
        solver_target_cartesian=np.asarray(new_arrival["evaluate"](boundary)).reshape(6).tolist(),
        terminal_command_clock_cone_rad=controls[-1].tolist(),
        solver_phase="handover",
        solver_status="NOT_RUN_EXPLICIT_CELL_HANDOVER",
        native_mission_accepted=False,
        handover_source_candidate_sha256=record["source_candidate_sha256"],
    )
    require(
        np.array_equal(np.asarray(result["controls_clock_cone_rad"]), controls),
        "Handover changed physical command angles",
    )
    host = validation.validate_candidate(
        result, ephemeris_path, stage="handover", mission_profile=profile
    )
    accepted = host.get("accepted_for_native_review") is True
    require(
        all(records.sha(path) == digest for path, digest in identities.items()),
        "Handover inputs changed",
    )
    record.update(
        status="HANDOVER_HOST_ACCEPTED" if accepted else "HANDOVER_REJECTED",
        reason="boundary_and_target_checked" if accepted else "retimed_boundary_host_check_failed",
        accepted_for_next_cell=accepted,
        source_host_validation=source_host,
        source_target_audit=source_target_audit,
        destination_target_audit=destination_target_audit,
        target_join=join,
        target_join_pass=join_pass,
        output_candidate_sha256=_digest(result),
        evidence_sha256=identities,
        wall_seconds=time.monotonic() - started,
    )
    return {"candidate": result, "host_validation": host, "record": record}
