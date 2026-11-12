"""ENC search, local refinement and native replay, in that order."""

import copy
import csv
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time
import uuid

from . import bridge, console, controls, ephemeris, objective, records, refinement
from . import report, optimiser, validation

ROOT = Path(__file__).resolve().parents[2]
SCHEMA = "long_fresh_workflow_v1"
DOMAIN = {"minimum_duration_days": 0.5, "maximum_duration_days": 600.0}
require, save, load = records.require, records.save, records.load


def profile(seed=1, evaluation_limit=200, mission_profile=None):
    require(type(seed) is int and 0 <= seed <= 2147483647, "Invalid random seed")
    require(
        type(evaluation_limit) is int
        and 40 <= evaluation_limit <= 1000000
        and (evaluation_limit % 2 == 0),
        "Use an even, complete-group cutoff >=40",
    )
    cells = [list(c) for c in ephemeris.campaign_cells(DOMAIN)]
    return {
        "schema": SCHEMA,
        "seed": seed,
        "evaluation_limit": evaluation_limit,
        "cold_start": True,
        "horizon_days": 600.0,
        "holds": 600,
        "objective": "minimum_time_at_declared_quality",
        "mission_profile": objective.validate_profile(
            mission_profile
            if mission_profile is not None
            else objective.reference_benchmark_profile()
        ),
        "native_score_contract": {
            "definition": objective.NATIVE_DEFINITION,
            "horizon_days": 600.0,
            "position_limit_m": 100000000.0,
            "velocity_limit_m_s": 100.0,
        },
        "native_score_role": "unchanged evaluation metric",
        "local_phases": ["restore", "quality"],
        "maximum_solver_calls": 1 + len(cells) * 3,
        "candidate_selection": "long_same_run_selection_v1",
        "arrival_cells_MJD": cells,
        "duration_bounds_days": [0.5, 600.0],
        "maximum_handover_retiming_days": 1e-06,
        "solver_worker_seconds": 1200.0,
        "solver_seconds_per_phase": 900.0,
        "solver_iterations_per_phase": 3000,
        "search_child_seconds": 720,
        "native_evaluation_seconds": 90,
        "native_global_position_km": 100000.0,
        "native_global_velocity_m_s": 100.0,
        "numerical_position_m": 100.0,
        "numerical_velocity_m_s": 0.001,
        "automatic_retry": False,
        "historical_controller_fallback": False,
        "continuation_domain": DOMAIN.copy(),
        "control_policy": controls.frozen_reference_policy(),
        "maximum_quality_solver_calls": len(cells) * 3,
        "same_cell_attempt_limit": 3,
        "optimisation_wall_budget_seconds": 3000.0,
        "declared_budget_progress_continuation": True,
    }


def files(directory):
    paths = list(Path(directory).rglob("*"))
    require(not any(path.is_symlink() for path in paths), "Linked evidence is not permitted")
    return sorted(path for path in paths if path.is_file())


def verify_seal(directory):
    binding = load(directory / "binding.json")
    outputs = [str(Path(path).relative_to(directory)) for path in binding["outputs"]]
    return records.verify(directory, binding["inputs"], outputs)


def implementation_files():
    return sorted((ROOT / "code/local_optimizer").glob("*.py")) + [
        ROOT / "run.py",
        ROOT / "code/local_optimizer/native_worker.ps1",
    ]


def prepare(
    installation_path,
    output,
    reference_commands,
    *,
    seed=1,
    evaluation_limit=200,
    mission_profile=None,
):
    """Prepare one new package; this does not run ENC or a local solver."""
    from . import bridge as native

    started = time.monotonic()
    reference = Path(reference_commands).resolve(strict=True)
    require(
        controls.build_policy(reference) == controls.frozen_reference_policy(),
        "Use the precise native-reference commands",
    )
    work = Path(output).absolute()
    require(not work.exists() and not work.is_symlink(), "Choose a new output directory")
    require(not work.resolve().is_relative_to(ROOT), "Keep raw runs outside source checkout")
    installation = load(installation_path)
    required = {
        "reference_inputs",
        "data_directory",
        "search_executable",
        "export_executable",
        "precision_executable",
        "search_sha256",
        "export_sha256",
        "precision_sha256",
        "plot_python",
    }
    require(set(installation) == required, "Installation must name exactly the required assets")
    for key in (
        "reference_inputs",
        "data_directory",
        "search_executable",
        "export_executable",
        "precision_executable",
        "plot_python",
    ):
        require(
            Path(installation[key]).is_absolute() and Path(installation[key]).exists(),
            "Missing absolute " + key,
        )
    config = profile(seed, evaluation_limit, mission_profile)
    local_budget = config.get(
        "optimisation_wall_budget_seconds",
        config["maximum_solver_calls"] * config["solver_worker_seconds"],
    )
    require(
        config["search_child_seconds"]
        + 3 * config["native_evaluation_seconds"]
        + local_budget
        + 900
        <= native.SESSION_SECONDS,
        "Worker session is shorter than the declared stage budgets and overhead",
    )
    from . import optimiser

    context = {
        key: config[key] for key in ("continuation_domain", "control_policy") if key in config
    }
    optimiser.protocol(
        config["mission_profile"], **context
    )  # Reject unsupported physics before any Windows job.
    config.update(run_id=uuid.uuid4().hex, installation=installation)
    work.mkdir(parents=True)
    (work / "_audit").mkdir()
    save(work / "workflow.json", config)
    try:
        native.configure(
            work / "native",
            installation["reference_inputs"],
            installation["search_executable"],
            installation["export_executable"],
            installation["precision_executable"],
            installation["data_directory"],
            installation["search_sha256"],
            installation["export_sha256"],
            installation["precision_sha256"],
            seed=seed,
            evaluation_limit=evaluation_limit,
            search_seconds=config["search_child_seconds"],
            **(
                {"terminal_policy": config["observer_terminal_policy"]}
                if "observer_terminal_policy" in config
                else {}
            ),
        )
        shutil.copy2(reference, work / "_audit/reference_commands.csv")
        inputs = implementation_files() + [work / "workflow.json", Path(installation_path)]
        inputs += files(ROOT / "case")
        inputs += [reference, work / "_audit/reference_commands.csv"]
        inputs += [p for p in (work / "native").rglob("*") if p.is_file()]
        save(work / "_audit/inputs.json", records.fingerprint(inputs))
        save(
            work / "preparation.json",
            {
                "status": "PREPARED_NOT_EXECUTED",
                "seconds": time.monotonic() - started,
                "new_enc_search": False,
                "local_solver_calls": 0,
            },
        )
    except Exception as error:
        from . import report

        elapsed = time.monotonic() - started
        save(work / "preparation_failure.json", {"reason": str(error), "seconds": elapsed})
        partial = {"status": "PREPARATION_FAILED", "accepted": False, "records": {}}
        save(
            work / "completion.json",
            {
                "status": "PREPARATION_FAILED",
                "scientific_accepted": False,
                "workflow_complete": False,
                "total_start_to_end_wall_seconds": elapsed,
                "six_metrics": report.metrics(partial, elapsed),
            },
        )
        raise
    return {
        "workspace": str(work),
        "launcher": str(work / "native/RUN_NATIVE_WORKER.cmd"),
        "status": "PREPARED_WINDOWS_WORKER_REQUIRED",
    }


def checked(work):
    from . import bridge as native

    config = load(work / "workflow.json")
    expected = profile(
        config.get("seed"), config.get("evaluation_limit"), config.get("mission_profile")
    )
    require(
        set(config) == set(expected) | {"run_id", "installation"}
        and all(config[k] == v and type(config[k]) is type(v) for k, v in expected.items()),
        "Workflow contract changed",
    )
    require(
        isinstance(config["run_id"], str) and len(config["run_id"]) == 32, "Run identity missing"
    )
    pins = load(work / "_audit/inputs.json")
    require(
        {str(p.resolve()) for p in implementation_files()} <= set(pins)
        and records.fingerprint(pins) == pins,
        "Prepared source/assets changed; preserve attempt",
    )
    native.verify(work / "native")
    require(
        controls.build_policy(Path(work) / "_audit/reference_commands.csv")
        == config["control_policy"],
        "Reference command limits changed",
    )
    return config


class Workflow:
    worker_entry = ROOT / "run.py"
    check_configuration = staticmethod(checked)

    def __init__(self, work):
        self.work = Path(work)
        self.console = console.Console(self.work)
        self.native_console = console.TextTail(self.work / "native/_transport/native_console.log")
        self.records = {}
        self.stage_times = {}
        self.selection = None

    def native_stage(self, name, schedule=None):
        from . import bridge as native

        start = time.monotonic()
        native.prepare_stage(self.work / "native", name, schedule=schedule)
        labels = {
            "cold_search": "[1/4] Fresh ENC search",
            "export": "[2/4] Export the saved controller: states and sail commands",
            "precise_baseline": "[3/4] Replay the starting schedule",
            "precise_candidate": "[4/4] Replay the optimised schedule from departure",
        }
        self.console.emit(labels[name])
        if name == "cold_search":
            self.console.emit("Waiting for the Windows worker. Keep both windows open.")
        waiting = time.monotonic()
        native.wait_ready(self.work / "native", name)
        self.stage_times[name + "_readiness_wait_seconds"] = time.monotonic() - waiting
        native.request(self.work / "native", name)

        def show_native():
            for line in self.native_console.read():
                self.console.emit(line, mirror=False)

        receipt = native.wait_result(self.work / "native", name, on_progress=show_native)
        show_native()
        self.console.emit(f"InTrance stage finished in {receipt['process_wall_seconds']:.2f} s.")
        stage = self.work / "native/stages" / name
        self.records[name + "_process"] = receipt
        self.stage_times[name + "_host_elapsed_seconds"] = time.monotonic() - start
        return stage

    def search(self):
        self.native_stage("cold_search")

    def export(self):
        from . import optimiser
        from .records import sha

        stage = self.native_stage("export")
        controller = stage / "fixture/mercury.eac"
        self.seed = optimiser.load_full_export(
            stage / "result",
            self.ephemeris,
            controller_sha256=sha(controller),
            run_id=self.config["run_id"],
            mission_profile=self.config["mission_profile"],
            **{
                key: self.config[key]
                for key in ("continuation_domain", "control_policy")
                if key in self.config
            },
        )
        out = self.work / "transfer"
        out.mkdir()
        save(out / "seed.json", self.seed)
        self.current = optimiser.baseline_candidate(
            self.seed, self.ephemeris, mission_profile=self.config["mission_profile"]
        )
        save(out / "candidate.json", self.current)
        self.schedule("baseline_schedule", self.current)

    def schedule(self, name, candidate):
        from . import validation

        output = self.work / name
        phase = "baseline" if name == "baseline_schedule" else self.selection["selected_phase"]
        # Only time-tagged commands and the original initial state are imported.
        validation.write_schedule(
            candidate, output, stage=phase, mission_profile=self.config["mission_profile"]
        )
        return output

    def audit(self, stage, candidate, label, *, phase):
        from . import validation

        output = self.work / (label + "_review")
        candidate_path = self.work / "_audit" / (label + "_candidate.json")
        save(candidate_path, candidate)
        binding = load(self.work / "_audit/inputs.json")
        binding.update(records.fingerprint(files(stage) + [candidate_path]))
        started = time.monotonic()
        info = validation.inspect_native_records(
            stage / "fixture",
            stage / "result",
            candidate,
            self.ephemeris,
            stage=phase,
            mission_profile=self.config["mission_profile"],
        )
        record = validation.numerical_review(info, output_directory=output)
        record["process_accepted"] = True
        record["process_receipt_accepted"] = True
        record["native_numerical_accepted"] = (
            record["recorded_native_numerical_checks_passed"] is True
        )
        record["native_mission_accepted"] = record["numerical_and_mission_gates_pass"] is True
        save(output / "native_comparison.json", record)
        records.seal(output, binding, [str(path.relative_to(output)) for path in files(output)])
        self.stage_times[label + "_validation_seconds"] = time.monotonic() - started
        self.records[label] = record
        numerical = record["native_numerical_accepted"]
        self.console.emit(
            f"{label.capitalize()} numerical comparison: {'passed' if numerical else 'failed'}."
        )
        encounter = record.get("encounter", {})
        if "position_miss_km" in encounter and "velocity_miss_m_s" in encounter:
            self.console.emit(
                f"Native arrival: {encounter['position_miss_km']:,.6f} km | "
                f"{encounter['velocity_miss_m_s']:.9f} m/s"
            )
        return info

    def baseline(self):
        stage = self.native_stage("precise_baseline", self.work / "baseline_schedule")
        info = self.audit(stage, self.current, "baseline", phase="baseline")
        if self.baseline_numerical_valid():
            from . import optimiser

            self.current = optimiser.with_precise_nodes(self.current, info)

    def baseline_numerical_valid(self):
        return self.records["baseline"].get("native_numerical_accepted") is True

    def replay(self):
        schedule = self.schedule("candidate_schedule", self.current)
        stage = self.native_stage("precise_candidate", schedule)
        self.audit(stage, self.current, "candidate", phase=self.selection["selected_phase"])

    def final_accepted(self):
        self.verify_completed()
        return self.records["candidate"].get("native_mission_accepted") is True

    def preflight(self):
        from . import ephemeris, validation

        self.config = self.check_configuration(self.work)
        self.console.setup(self.config)
        self.ephemeris = Path(self.config["installation"]["data_directory"]) / "bin1600_600.405"
        descriptors = ephemeris.descriptors(self.ephemeris, self.config["continuation_domain"])
        require(
            [d["cell_MJD"] for d in descriptors] == self.config["arrival_cells_MJD"],
            "Long domain metadata changed",
        )
        # Later cells are checked when entered, not all symbolically built up front.
        review = validation.validate_target_model(
            self.ephemeris,
            ephemeris.build_cell(self.ephemeris, descriptors[0]),
            continuation_domain=self.config["continuation_domain"],
        )
        save(
            self.work / "_audit/target_cells.json",
            {
                "declared_cells": self.config["arrival_cells_MJD"],
                "initial_cell_review": review,
                "later_cells_checked_on_entry": True,
            },
        )
        self.checkpoint_count = 0
        self.work.joinpath("checkpoints").mkdir()

    def checkpoint(self, incumbent, selected_key, ledger, reason, *, native_accepted=False):
        self.checkpoint_count += 1
        record = {
            "schema": "long_checkpoint_v1",
            "sequence": self.checkpoint_count,
            "run_id": self.config["run_id"],
            "selected_stage_key": selected_key,
            "candidate": incumbent["candidate"],
            "host_validation": incumbent["host_validation"],
            "solver_result": incumbent["result"],
            "attempts": ledger,
            "reason": reason,
            "remaining_optimisation_seconds": self.budget.remaining(),
            "native_accepted": native_accepted,
            "automatic_resume_authorized": False,
        }
        directory = self.work / "checkpoints"
        name = f"checkpoint_{self.checkpoint_count:04d}.json"
        temp = directory / (name + ".tmp")
        with temp.open("x", encoding="utf-8") as stream:
            json.dump(record, stream, indent=2, allow_nan=False)
            stream.write("\n")
            stream.flush()
            os.fsync(stream.fileno())
        require(not (directory / name).exists(), "Checkpoint collision")
        os.replace(temp, directory / name)
        pointer = directory / ("latest_" + str(self.checkpoint_count) + ".tmp")
        save(pointer, {"checkpoint": name, "sha256": records.sha(directory / name)})
        os.replace(pointer, directory / "latest.json")

    def invoke(self, phase, index, attempt, candidate, allocation, origin):
        key = "restore" if phase == "restore" else refinement.key(index, attempt)
        label = (
            "Restore arrival quality"
            if phase == "restore"
            else f"Reduce flight time: interval {index + 1}, attempt {attempt + 1}"
        )
        bounds = candidate["duration_bounds_days"]
        self.console.emit(
            f"{label} | starting TOF {candidate['duration_days']:.6f} days | "
            f"interval {bounds[0]:g}..{bounds[1]:g} days"
        )
        self.console.emit(
            f"Solver limit {allocation['solver_seconds']:g} s | "
            f"local budget remaining {self.budget.remaining():.0f} s. "
            "Iteration values are trials; continuous arrival checks follow the solve."
        )
        directory = self.work / ("local_" + key)
        directory.mkdir()
        save(directory / "input_candidate.json", candidate)
        context = {
            "phase": phase,
            "cell_index": index,
            "attempt_index": attempt,
            "allocation": allocation,
            "origin": origin,
            "stage_key": key,
        }
        save(directory / "context.json", context)
        inputs = load(self.work / "_audit/inputs.json")
        origin_files = [self.work / "baseline_review/native_comparison.json"]
        if phase == "restore":
            origin_files += [
                self.work / "transfer/candidate.json",
                self.work / "native/stages/precise_baseline/result/bridge_states.csv",
            ]
        if origin["kind"] != "baseline":
            previous = self.work / ("local_" + origin["previous_stage"])
            verify_seal(previous)
            origin_files += [previous / "result.json", previous / "binding.json"]
        if origin["kind"] == "handover":
            moved = self.work / origin["handover_directory"]
            verify_seal(moved)
            origin_files += [
                moved / n for n in ("input_candidate.json", "result.json", "binding.json")
            ]
        inputs.update(
            records.fingerprint(
                origin_files + [directory / "input_candidate.json", directory / "context.json"]
            )
        )
        save(directory / "inputs.json", inputs)
        env = dict(
            os.environ,
            OMP_NUM_THREADS="1",
            OPENBLAS_NUM_THREADS="1",
            MKL_NUM_THREADS="1",
            PYTHONNOUSERSITE="1",
            PYTHONUNBUFFERED="1",
        )
        started = time.monotonic()
        timeout = False
        with (
            (directory / "stdout.txt").open("xb") as out,
            (directory / "stderr.txt").open("xb") as err,
        ):
            remaining = self.budget.remaining()
            if remaining < 1.0:
                save(
                    directory / "process.json",
                    {
                        "return_code": None,
                        "timed_out": False,
                        "not_started": True,
                        "owned_child_cleanup_verified": True,
                        "wall_seconds": time.monotonic() - started,
                    },
                )
                records.seal(
                    directory,
                    inputs,
                    [
                        "input_candidate.json",
                        "context.json",
                        "process.json",
                        "stdout.txt",
                        "stderr.txt",
                    ],
                )
                raise refinement.BudgetExpired(
                    {
                        "stage_key": key,
                        "solver_status": "NOT_STARTED_BUDGET",
                        "solver_calls": 0,
                        "trial_admitted": False,
                    }
                )
            dispatch_timeout = min(allocation["worker_seconds"], remaining)
            child = subprocess.Popen(
                [
                    sys.executable,
                    "-B",
                    str(self.worker_entry),
                    "solve-worker",
                    str(self.work),
                    "--stage",
                    key,
                ],
                cwd=directory,
                env=env,
                stdout=out,
                stderr=err,
            )
            try:
                code = console.wait_for_solver(
                    child,
                    dispatch_timeout,
                    self.console,
                    directory / "stdout.txt",
                    phase=phase,
                    time_scale_days=objective.dynamical_time_days(self.config["mission_profile"]),
                    solver_limit=allocation["solver_seconds"],
                )
            except subprocess.TimeoutExpired:
                timeout = True
            finally:
                if child.poll() is None:
                    child.kill()
                code = child.wait(timeout=10)
                save(
                    directory / "process.json",
                    {
                        "return_code": code,
                        "timed_out": timeout,
                        "owned_child_cleanup_verified": child.poll() is not None,
                        "remaining_at_dispatch_seconds": remaining,
                        "dispatch_timeout_seconds": dispatch_timeout,
                        "wall_seconds": time.monotonic() - started,
                    },
                )
        outputs = [p.name for p in directory.iterdir() if p.is_file() and p.name != "inputs.json"]
        records.seal(directory, inputs, outputs)
        if timeout:
            # No incomplete vector is promoted. The chain may retain its previous point.
            raise refinement.BudgetExpired(
                {
                    "stage_key": key,
                    "solver_status": "HOST_WORKER_DEADLINE",
                    "stop_reason": (
                        "optimisation_wall_budget_exhausted"
                        if self.budget.remaining() < 1.0
                        else "local_worker_deadline_retained_checkpoint"
                    ),
                    "solver_calls": int((directory / "solver_call.json").exists()),
                    "trial_admitted": False,
                }
            )
        require(code == 0, "Long local worker failed; evidence preserved")
        result = load(directory / "result.json")
        self.records[key] = result
        self.console.local_result(key, result)
        return result

    def optimise(self):
        self.budget = refinement.Budget(self.config["optimisation_wall_budget_seconds"])
        try:
            self._optimise_with_budget()
        finally:
            self.budget.stop()
            measured = {
                "wall_seconds": self.budget.elapsed(),
                "budget_seconds": self.budget.limit,
                "overrun_seconds": max(0.0, self.budget.elapsed() - self.budget.limit),
                "scope": (
                    "restoration, refinement, handovers and checkpoints; final "
                    "native replay excluded"
                ),
                "overrun_not_hidden": True,
            }
            save(self.work / "optimisation_budget.json", measured)
            self.stage_times["optimisation_phase_wall_seconds"] = measured["wall_seconds"]

    def _optimise_with_budget(self):
        restored = self.invoke(
            "restore", 0, 0, self.current, self.budget.allocation(), {"kind": "baseline"}
        )
        require(refinement.admitted(restored), "No same-run guarded restoration")
        origins = {0: {"kind": "restoration", "previous_stage": "restore"}}

        def solve(index, attempt, candidate, allocation):
            origin = (
                origins[index]
                if attempt == 0
                else {
                    "kind": "same_cell_retry",
                    "previous_stage": refinement.key(index, attempt - 1),
                }
            )
            return self.invoke("quality", index, attempt, candidate, allocation, origin)

        def handover(index, candidate, cell, previous_key):
            self.console.emit(
                f"Reached {candidate['duration_days']:.6f} days at an interval boundary. "
                "Checking the move into the next earlier interval."
            )
            directory = self.work / f"handover_{index + 1:03d}"
            directory.mkdir()
            require(
                candidate == self.records[previous_key]["candidate"], "Handover predecessor differs"
            )
            inputs = load(self.work / "_audit/inputs.json")
            previous = self.work / ("local_" + previous_key)
            inputs.update(
                records.fingerprint([previous / "result.json", previous / "binding.json"])
            )
            save(directory / "input_candidate.json", candidate)
            save(
                directory / "context.json",
                {"cell_index": index, "next_cell_MJD": cell, "previous_stage": previous_key},
            )
            inputs.update(
                records.fingerprint(
                    [directory / "input_candidate.json", directory / "context.json"]
                )
            )
            save(directory / "inputs.json", inputs)
            start = time.monotonic()
            timeout = False
            with (
                (directory / "stdout.txt").open("xb") as out,
                (directory / "stderr.txt").open("xb") as err,
            ):
                remaining = self.budget.remaining()
                if remaining <= 0.0:
                    result = {
                        "candidate": None,
                        "host_validation": None,
                        "record": {
                            "status": "HANDOVER_NOT_STARTED",
                            "reason": "optimisation_wall_budget_exhausted",
                        },
                    }
                    save(directory / "deadline_result.json", result)
                    save(
                        directory / "process.json",
                        {
                            "not_started": True,
                            "return_code": None,
                            "timed_out": False,
                            "owned_child_cleanup_verified": True,
                            "wall_seconds": time.monotonic() - start,
                        },
                    )
                    records.seal(
                        directory,
                        inputs,
                        [
                            p.name
                            for p in directory.iterdir()
                            if p.is_file() and p.name != "inputs.json"
                        ],
                    )
                    return result
                deadline = min(120.0, remaining)
                child = subprocess.Popen(
                    [
                        sys.executable,
                        "-B",
                        str(self.worker_entry),
                        "handover-worker",
                        str(self.work),
                        "--directory",
                        directory.name,
                    ],
                    cwd=directory,
                    stdout=out,
                    stderr=err,
                    env=dict(
                        os.environ,
                        OMP_NUM_THREADS="1",
                        OPENBLAS_NUM_THREADS="1",
                        MKL_NUM_THREADS="1",
                        PYTHONNOUSERSITE="1",
                    ),
                )
                try:
                    code = child.wait(timeout=deadline)
                except subprocess.TimeoutExpired:
                    timeout = True
                finally:
                    if child.poll() is None:
                        child.kill()
                    code = child.wait(timeout=10)
                    save(
                        directory / "process.json",
                        {
                            "return_code": code,
                            "timed_out": timeout,
                            "owned_child_cleanup_verified": child.poll() is not None,
                            "wall_seconds": time.monotonic() - start,
                            "dispatch_timeout_seconds": deadline,
                        },
                    )
            if timeout:
                result = {
                    "candidate": None,
                    "host_validation": None,
                    "record": {
                        "status": "HANDOVER_DEADLINE",
                        "reason": "handover_deadline_retained_checkpoint",
                    },
                }
                save(directory / "deadline_result.json", result)
            else:
                require(code == 0, "Handover worker failed; checkpoint retained but campaign stops")
                result = load(directory / "result.json")
            records.seal(
                directory,
                inputs,
                [p.name for p in directory.iterdir() if p.is_file() and p.name != "inputs.json"],
            )
            origins[index] = {
                "kind": "handover",
                "previous_stage": previous_key,
                "handover_directory": directory.name,
            }
            return result

        self.current, self.selection = refinement.run_with_budget(
            restored,
            self.config["arrival_cells_MJD"],
            solve,
            handover,
            self.checkpoint,
            self.budget,
        )
        save(self.work / "selection.json", self.selection)

    def verify_completed(self):
        from . import bridge

        self.check_configuration(self.work)
        require(
            load(self.work / "selection.json") == self.selection
            and self.current == self.records[self.selection["selected_stage_key"]]["candidate"],
            "Selected point changed",
        )
        for stage in bridge.STAGES:
            bridge.inspect(self.work / "native/stages" / stage)
        quality = list(self.work.glob("local_quality_*"))
        require(
            len(quality)
            == self.selection["quality_attempts"]
            <= self.config["maximum_quality_solver_calls"],
            "Unrecorded local attempts",
        )
        for directory in (
            [
                self.work / "baseline_review",
                self.work / "candidate_review",
                self.work / "local_restore",
            ]
            + quality
            + list(self.work.glob("handover_*"))
        ):
            verify_seal(directory)
        pointer = load(self.work / "checkpoints/latest.json")
        require(
            records.sha(self.work / "checkpoints" / pointer["checkpoint"]) == pointer["sha256"],
            "Checkpoint changed",
        )


def execute(backend):
    backend.preflight()
    backend.search()
    backend.export()
    backend.baseline()
    require(backend.baseline_numerical_valid(), "Baseline numerical transfer rejected")
    backend.optimise()
    backend.replay()
    accepted = backend.final_accepted()
    incumbent = backend.records[backend.selection["selected_stage_key"]]
    backend.checkpoint(
        incumbent,
        backend.selection["selected_stage_key"],
        backend.selection["cells"],
        "final_native_accepted" if accepted else "final_native_rejected",
        native_accepted=accepted,
    )
    return accepted


def checked_stage(stage):
    if stage == "restore":
        return stage
    parts = stage.split("_")
    require(
        len(parts) == 3 and parts[0] == "quality" and parts[1].isdigit() and parts[2].isdigit(),
        "Invalid stage key",
    )
    require(stage == refinement.key(int(parts[1]) - 1, int(parts[2]) - 1), "Noncanonical stage key")
    return stage


def handover_worker(work, directory_name):
    from . import refinement

    work = Path(work).resolve()
    config = checked(work)
    require(
        directory_name.startswith("handover_") and directory_name[9:].isdigit(),
        "Invalid handover directory",
    )
    directory = work / directory_name
    require(not (directory / "result.json").exists(), "Handover already attempted")
    context = load(directory / "context.json")
    index = context["cell_index"]
    require(
        type(index) is int
        and 0 < index < len(config["arrival_cells_MJD"])
        and directory_name == f"handover_{index + 1:03d}"
        and context["next_cell_MJD"] == config["arrival_cells_MJD"][index],
        "Handover domain differs",
    )
    inputs = load(directory / "inputs.json")
    require(records.fingerprint(inputs) == inputs, "Handover inputs changed")
    previous_dir = work / ("local_" + checked_stage(context["previous_stage"]))
    verify_seal(previous_dir)
    previous = load(previous_dir / "result.json")
    candidate = load(directory / "input_candidate.json")
    require(
        refinement.admitted(previous)
        and candidate == previous["candidate"]
        and candidate["arrival_cell_MJD"] == config["arrival_cells_MJD"][index - 1]
        and candidate["run_id"] == config["run_id"]
        and candidate["control_policy"] == config["control_policy"]
        and candidate["continuation_domain"] == config["continuation_domain"],
        "Handover predecessor differs",
    )
    result = refinement.prepare_handover(
        candidate,
        Path(config["installation"]["data_directory"]) / "bin1600_600.405",
        context["next_cell_MJD"],
    )
    save(directory / "result.json", result)
    require(records.fingerprint(inputs) == inputs, "Handover inputs changed during review")


def solve_worker(work, stage):
    from . import optimiser, bridge

    work = Path(work).resolve()
    config = checked(work)
    checked_stage(stage)
    directory = work / ("local_" + stage)
    context = load(directory / "context.json")
    phase, index, attempt = context["phase"], context["cell_index"], context["attempt_index"]
    expected = "restore" if phase == "restore" else refinement.key(index, attempt)
    require(
        stage == expected
        and phase in ("restore", "quality")
        and 0 <= index < len(config["arrival_cells_MJD"]),
        "Worker stage differs",
    )
    require(not (directory / "attempt.json").exists(), "Solver attempt already used")
    inputs = load(directory / "inputs.json")
    require(records.fingerprint(inputs) == inputs, "Worker inputs changed")
    candidate = load(directory / "input_candidate.json")
    require(
        candidate["run_id"] == config["run_id"]
        and candidate["mission_profile"] == config["mission_profile"]
        and candidate["control_policy"] == config["control_policy"]
        and candidate["continuation_domain"] == config["continuation_domain"],
        "Candidate contract differs",
    )
    for stage_id in ("cold_search", "export", "precise_baseline"):
        bridge.inspect(work / "native/stages" / stage_id)
    require(
        candidate["origin_controller_sha256"]
        == records.sha(work / "native/stages/export/fixture/mercury.eac"),
        "Controller lineage differs",
    )
    origin = context["origin"]
    verify_seal(work / "baseline_review")
    baseline = load(work / "baseline_review/native_comparison.json")
    require(
        baseline.get("native_numerical_accepted") is True
        and baseline.get("process_receipt_accepted") is True,
        "Restoration needs the same-run audited numerical baseline",
    )
    if phase == "restore":
        from . import transfer

        require(
            origin == {"kind": "baseline"} and index == attempt == 0,
            "Unexpected restoration origin",
        )
        expected_candidate = load(work / "transfer/candidate.json")
        with (work / "native/stages/precise_baseline/result/bridge_states.csv").open(
            newline=""
        ) as stream:
            rows = list(csv.DictReader(stream))
        expected_candidate["nodes_cartesian"] = [
            [float(row[name]) for name in transfer.CARTESIAN_FIELDS] for row in rows
        ]
        expected_candidate["starting_nodes"] = "same-run precise native baseline"
        require(
            candidate == expected_candidate, "Restoration input is not the precise fresh baseline"
        )
    if phase == "quality":
        previous_dir = work / ("local_" + checked_stage(origin["previous_stage"]))
        verify_seal(previous_dir)
        previous = load(previous_dir / "result.json")
        require(refinement.admitted(previous), "Prior point is not admitted")
        if origin["kind"] == "handover":
            moved_dir = work / origin["handover_directory"]
            verify_seal(moved_dir)
            moved = load(moved_dir / "result.json")
            require(
                load(moved_dir / "input_candidate.json") == previous["candidate"]
                and candidate == moved["candidate"]
                and moved["host_validation"]["accepted_for_native_review"] is True,
                "Handover identity differs",
            )
        else:
            require(
                origin["kind"] in ("restoration", "same_cell_retry")
                and candidate == previous["candidate"],
                "Restart input differs",
            )
        require(
            candidate["arrival_cell_MJD"] == config["arrival_cells_MJD"][index],
            "Input cell differs",
        )
    allocation = context["allocation"]
    require(
        60.0 <= allocation["solver_seconds"] <= 900.0
        and allocation["solver_seconds"] < allocation["worker_seconds"] <= 1200.0,
        "Invalid worker allocation",
    )
    save(
        directory / "attempt.json", {"phase": phase, "stage_key": stage, "maximum_solver_calls": 1}
    )
    result = optimiser.solve_candidate(
        candidate,
        Path(config["installation"]["data_directory"]) / "bin1600_600.405",
        phase=phase,
        before_solver_call=lambda record: save(directory / "solver_call.json", record),
        mission_profile=config["mission_profile"],
        solver_wall_budget_seconds=allocation["solver_seconds"],
    )
    save(directory / "result.json", result)
    require(records.fingerprint(inputs) == inputs, "Inputs changed during solve")


def run(workspace, *, started=None):
    """One monotonic outer timer includes waiting, all stages, reports and cleanup."""
    from . import bridge as native, report

    includes_preparation = started is not None
    started = time.monotonic() if started is None else started
    work = Path(workspace).resolve()
    require(not (work / "attempt.json").exists(), "Attempt already used; choose a new workspace")
    save(
        work / "attempt.json",
        {
            "started_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            "timing_basis": (
                "host monotonic, includes readiness waits; final timing " "receipt write excluded"
            ),
        },
    )
    backend = Workflow(work)
    status, error, accepted = "INCOMPLETE", None, False
    try:
        accepted = execute(backend)
        status = "FRESH_WORKFLOW_NATIVE_ACCEPTED" if accepted else "NATIVE_REJECTED"
    except (Exception, KeyboardInterrupt) as exc:
        status, error = "STOPPED_NO_RETRY", type(exc).__name__ + ": " + str(exc)
    try:
        cleanup = native.stop_and_wait(
            work / "native",
            request_stop=status != "FRESH_WORKFLOW_NATIVE_ACCEPTED",
            reason="host_admission_stop" if status == "STOPPED_NO_RETRY" else "host_finalization",
        )
    except Exception as exc:
        cleanup = {
            "terminal_acknowledged": False,
            "owned_child_cleanup_verified": False,
            "reason": str(exc),
        }
    if (
        cleanup.get("terminal_acknowledged") is not True
        or cleanup.get("owned_child_cleanup_verified") is not True
    ):
        accepted, status = False, "CLEANUP_UNVERIFIED"
    result = {
        "schema": SCHEMA,
        "status": status,
        "accepted": accepted,
        "error": error,
        "records": backend.records,
        "stage_times": backend.stage_times,
        "cleanup": cleanup,
        "mission_profile": getattr(backend, "config", {}).get(
            "mission_profile", objective.default_profile()
        ),
        "selection": getattr(backend, "selection", None),
        "transport_failure": cleanup.get("worker", {}).get("state") == "FAILED",
        "error_origin": (
            "native_worker"
            if cleanup.get("worker", {}).get("state") == "FAILED"
            else "host_solver_or_validation" if error else None
        ),
        "local_processes": {
            directory.name.removeprefix("local_"): load(directory / "process.json")
            for directory in sorted(work.glob("local_*"))
            if (directory / "process.json").exists()
        },
        "fresh_enc_requested": True,
        "global_optimum_claim": False,
        "equal_budget_superiority_claim": False,
        "original_files_modified": False,
    }
    save(work / "scientific_result.json", result)
    report_error = None
    try:
        report.generate(work, result)
        if accepted:
            backend.verify_completed()
        report_files = files(work / "report")
        save(
            work / "_audit/report_binding.json",
            {
                "inputs": records.fingerprint([work / "scientific_result.json"]),
                "outputs": records.fingerprint(report_files),
            },
        )
    except Exception as exc:
        report_error = type(exc).__name__ + ": " + str(exc)
    try:
        six_metrics = report.metrics(result)
    except Exception as exc:
        report_error = (
            (report_error or "") + " Metrics failed: " + type(exc).__name__ + ": " + str(exc)
        )
        six_metrics = {
            "final_position_miss_km": None,
            "final_velocity_miss_m_s": None,
            "final_flight_time_days": None,
            "final_native_fitness": None,
            "final_own_fitness": None,
            "total_start_to_end_wall_seconds": None,
            "metrics_error": str(exc),
        }
    elapsed = time.monotonic() - started
    six_metrics["total_start_to_end_wall_seconds"] = elapsed
    completion = {
        "status": status if report_error is None else "REPORT_INCOMPLETE",
        "scientific_status": status,
        "scientific_accepted": accepted,
        "workflow_complete": bool(accepted and report_error is None),
        "report_complete": report_error is None,
        "report_error": report_error,
        "total_start_to_end_wall_seconds": elapsed,
        "preparation_included_in_total": includes_preparation,
        "separate_preparation_seconds": (
            load(work / "long_preparation.json")["seconds"]
            if (work / "long_preparation.json").exists()
            else (
                load(work / "preparation.json")["seconds"]
                if (work / "preparation.json").exists()
                else None
            )
        ),
        "total_includes": [
            "preflight",
            "fresh ENC",
            "adapter",
            "all local calls",
            "native replay",
            "independent checks",
            "readiness waits",
            "report and plot attempt",
            "owned cleanup",
        ],
        "excluded": ["one-time installation", "final timing receipt assembly/write"]
        + (
            []
            if includes_preparation
            else ["earlier package preparation; see separate_preparation_seconds"]
        ),
        "nested_stage_times_must_not_be_added": True,
        "six_metrics": six_metrics,
    }
    save(work / "completion.json", completion)
    console.finish(
        backend.console if hasattr(backend, "console") else console.Console(work), completion
    )
    return completion
