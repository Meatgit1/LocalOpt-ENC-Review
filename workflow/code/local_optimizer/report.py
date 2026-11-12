"""Report six experiment metrics and plot the saved InTrance results."""

import csv
import json
import math
import re
from pathlib import Path
import subprocess
from . import records, objective


def _row(position_km, velocity_m_s, duration, basis, profile, *, gates=None):
    own = objective.evaluate(position_km * 1000.0, velocity_m_s, duration, profile, **gates or {})
    return {
        "schema": "mission_scored_endpoint_v2",
        "evaluation_basis": basis,
        "metrics": {
            "arrival_position_km": position_km,
            "arrival_position_m": position_km * 1000.0,
            "arrival_velocity_m_s": velocity_m_s,
            "flight_time_days": duration,
        },
        "native_fitness": objective.native_endpoint_fitness(
            position_km * 1000.0, velocity_m_s, duration
        ),
        "own_fitness": own["own_fitness"],
        "own_evaluation": own,
        "scores_establish_acceptance": False,
    }


def metrics(science, elapsed=None):
    rows = []
    profile = objective.validate_profile(science.get("mission_profile"))
    search = science["records"].get("cold_search_process", {}).get("search_summary")
    if search and search.get("selected_endpoint"):
        selected = search["selected_endpoint"]
        row = _row(
            selected["position_miss_km"],
            selected["velocity_miss_m_s"],
            selected["flight_time_days"],
            "selected_native_encounter",
            profile,
        )
        row.update(
            label="native_selected",
            printed_precision_only=True,
            reported_native_fitness=search["native_fitness"],
            native_mission_accepted=None,
            whole_workflow_accepted=science["accepted"],
        )
        rows.append(row)
    for key in ("baseline", "candidate"):
        record = science["records"].get(key)
        if not record:
            continue
        encounter = record["encounter"]
        row = _row(
            encounter["position_miss_km"],
            encounter["velocity_miss_m_s"],
            record["duration_days"],
            "fixed_endpoint",
            profile,
            gates={
                "native_validated": record.get("process_receipt_accepted") is True,
                "numerical_validated": record.get("native_numerical_accepted") is True,
                "physical_validated": record.get("scientific_records_bound") is True
                and record.get("control_regularity_pass", True) is True,
            },
        )
        require_match = row["own_evaluation"]["fully_validated_quality_eligible"]
        if (record.get("native_mission_accepted") is True) != require_match:
            raise ValueError("Native acceptance and declared local-quality gates disagree")
        row.update(
            label=key,
            native_mission_accepted=record.get("native_mission_accepted") is True,
            whole_workflow_accepted=science["accepted"],
        )
        rows.append(row)
    selected = next((row for row in rows if row["label"] == "candidate"), None)
    baseline = next((row for row in rows if row["label"] == "baseline"), None)
    return {
        "schema": "fresh_six_metrics_v2",
        "status": science["status"],
        "mission_profile": profile,
        "selection": science.get("selection"),
        "final_position_miss_km": (
            None if selected is None else selected["metrics"]["arrival_position_km"]
        ),
        "final_velocity_miss_m_s": (
            None if selected is None else selected["metrics"]["arrival_velocity_m_s"]
        ),
        "final_flight_time_days": (
            None if selected is None else selected["metrics"]["flight_time_days"]
        ),
        "total_start_to_end_wall_seconds": elapsed,
        "final_native_fitness": None if selected is None else selected["native_fitness"],
        "final_own_fitness": None if selected is None else selected["own_fitness"],
        "reference_benchmark_time_improved": (
            None
            if selected is None
            else selected["own_evaluation"].get("reference_benchmark_time_improved")
        ),
        "trajectories": rows,
        "native_search_selected_encounter": science["records"]
        .get("cold_search_process", {})
        .get("search_summary"),
        "own_score_increased_vs_full_horizon": (
            None
            if selected is None
            or baseline is None
            or selected["own_fitness"] is None
            or (baseline["own_fitness"] is None)
            else selected["own_fitness"] > baseline["own_fitness"] + 1e-08
        ),
        "both_histories_pass_arrival_gates": (
            None
            if selected is None or baseline is None
            else selected["native_mission_accepted"] and baseline["native_mission_accepted"]
        ),
        "missing_final_metrics_are_not_zero": True,
        "baseline_is_full_horizon_not_selected_native_encounter": True,
        "own_score_does_not_measure_control_slew_or_uncertainty": True,
    }


def _budget_scope(work, science, config):
    """Describe the recorded allocation, not a hardcoded historical profile."""
    configured = config.get("optimisation_wall_budget_seconds")
    allocation = (science.get("selection") or {}).get("optimisation_wall_budget_seconds")
    budget_path = work / "optimisation_budget.json"
    if allocation is None and budget_path.exists():
        allocation = json.loads(budget_path.read_text())["budget_seconds"]
    for value in (configured, allocation):
        if value is not None and (
            type(value) not in (int, float) or not math.isfinite(value) or value <= 0
        ):
            raise ValueError("Invalid recorded local-phase budget")
    if configured is not None and allocation is not None and (allocation > configured):
        raise ValueError("Recorded local allocation exceeds its configured maximum")
    maximum = (
        f"Local budget: up to {configured:g} s. "
        if configured is not None
        else "Local budget: not recorded. "
    )
    admitted = (
        f"Available for this local stage: {allocation:g} s. "
        if allocation is not None
        else "Available local time: not recorded. "
    )
    return (
        "One seed and evaluation cutoff. "
        + maximum
        + admitted
        + (
            "This covers model building, restoration, time solves, handovers, checks "
            "and saving. Each call uses the remaining allowance. Convergence or a "
            "failed check can end it early."
        )
    )


def generate(workspace, science):
    work = Path(workspace)
    output = work / "report"
    output.mkdir()
    measured = metrics(science)
    config_path = work / "workflow.json"
    config = json.loads(config_path.read_text()) if config_path.exists() else {}
    records.save(output / "metrics_before_timing.json", measured)
    lines = [
        "# Earth-Mercury run",
        "",
        "Status: " + science["status"],
        "",
        "Fresh ENC -> exported flight -> local steering/time changes -> InTrance replay.",
        "",
        (
            "| Recorded trajectory | Flight [days] | Position miss [km] |"
            " Velocity miss [m/s] | Native score | Time score | Arrival + "
            "numerical checks |"
        ),
        "|---|---:|---:|---:|---:|---:|---|",
    ]
    for row in measured["trajectories"]:
        raw = row["metrics"]
        native_score = row.get("reported_native_fitness", row["native_fitness"])
        gates = (
            "not independently checked"
            if row["native_mission_accepted"] is None
            else "pass" if row["native_mission_accepted"] else "not passed"
        )
        own_score = "pending" if row["own_fitness"] is None else f"{row['own_fitness']:.6f}"
        label = {
            "native_selected": "ENC selected point",
            "baseline": "Starting 600-day flight",
            "candidate": "Local schedule (native replay)",
        }[row["label"]]
        lines.append(
            f"| {label} | {raw['flight_time_days']:.6f} | "
            f"{raw['arrival_position_km']:.3f} | {raw['arrival_velocity_m_s']:.6f} | "
            f"{native_score:.6f} | {own_score} | {gates} |"
        )
    if not measured["trajectories"]:
        lines += [
            "",
            "No checked native trajectory was returned. Missing values stay unavailable.",
        ]
    lines += [
        "",
        (
            "The starting flight is the full 600-day ENC history. ENC selected point is its "
            "encounter."
        ),
        (
            "The ENC selected point uses the normal CSV and search score. Its time score "
            "stays pending. The other rows use precise endpoint records."
        ),
        (
            "The two scores have different definitions. Time is ranked after the "
            "recorded arrival, numerical and control checks."
        ),
        f"Local arrival limits: {measured['mission_profile']['position_target_m']:g} m / "
        f"{measured['mission_profile']['velocity_target_m_s']:g} m/s. "
        "Global ENC limits: 100,000 km / 100 m/s. "
        "Numerical-agreement limits: 100 m / 0.001 m/s.",
        ("Path checks use the integration samples. Navigation errors are outside " "this model."),
        "",
        "## Computer time",
        "",
        (
            "[completion.json](../completion.json) has the final numbers and full "
            "elapsed time, including waits, reporting and cleanup. Solver time is "
            "already inside that total."
        ),
        ("Failed attempts keep their own partial results and elapsed time."),
        "",
        "## Local budget",
        "",
        _budget_scope(work, science, config),
        "Native-checked time improvement over 576 days (None = unavailable): "
        f"{measured['reference_benchmark_time_improved']}.",
    ]
    if science.get("error"):
        lines += [
            "",
            "Stopped: " + science["error"],
            "Error from: " + str(science.get("error_origin", "not classified")) + ".",
        ]
    if science.get("selection"):
        selection = science["selection"]
        lines += [
            "",
            "## What was kept?",
            "",
            f"Selected phase: **{selection['selected_phase']}** ({selection['reason']}).",
            f"Solver status: {selection['refinement_status']}; "
            f"host checks passed: {selection['refinement_admitted']}.",
            "Solver converged: "
            f"{selection.get('selected_solver_converged', 'not recorded')}; "
            f"kept as: {selection.get('selected_admission', 'legacy')}.",
            (
                "A time-limited solve can still return a usable point. Its convergence "
                "status is shown above."
            ),
            (
                "Keeping restoration means no shorter checked candidate was selected. The "
                "native outcome is in the result table."
            ),
        ]
        if selection.get("cells"):
            lines += [
                "",
                (
                    "| Time interval attempt | Selected flight [days] | Solver "
                    "status | Host checks passed |"
                ),
                "|---|---:|---|---|",
            ]
            for item in selection["cells"]:
                duration = (
                    f"{item['selected_duration_days']:.6f}"
                    if "selected_duration_days" in item
                    else "retained checkpoint"
                )
                lines.append(
                    f"| {item['cell_index'] + 1} | {duration} | {item['solver_status']} | "
                    f"{item['trial_admitted']} |"
                )
        lines += [
            "",
            "Consecutive sail-normal commands use the reference-derived limit.",
            (
                "Commands change at hold boundaries. Sail turning dynamics are outside this"
                " model."
            ),
            "checkpoints/ keeps the admitted local candidates and the final native outcome.",
        ]
    lines += [
        "",
        "## Plots",
        "",
        "trajectory_quality: native paths and the two scores.",
        "controls: held clock/cone and changes between neighbouring command normals.",
        "refinement_progress: kept times, attempted cells and the last available "
        "Ipopt log. Time solves show approximate TOF; restoration shows arrival "
        "cost. The rows are reported iterates, and inf_pr uses the submitted NLP "
        "units.",
        "numerical_checks: native-RK4-256, RK4-128-RK4-256 and nodes-RK4-256 "
        "differences in m and m/s, with separate limits. 128/256 mean RK4 substeps "
        "per command. Missing histories use labelled saved maxima where available.",
        "See [plot_index.json](plot_index.json) for the figure files, source "
        "records and any missing part.",
    ]
    (output / "README.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    if any(row["label"] != "native_selected" for row in measured["trajectories"]) or any(
        work.glob("local_*/result.json")
    ):
        config = json.loads((work / "workflow.json").read_text())
        renderer = Path(__file__).resolve().parents[2] / "run.py"
        with (
            (output / "plot.stdout.txt").open("xb") as stdout,
            (output / "plot.stderr.txt").open("xb") as stderr,
        ):
            completed = subprocess.run(
                [config["installation"]["plot_python"], str(renderer), "plot", str(work)],
                stdout=stdout,
                stderr=stderr,
                timeout=60,
                check=False,
            )
        if completed.returncode != 0:
            raise RuntimeError("Plot failed; see preserved plot.stderr.txt")
        index_path = output / "plot_index.json"
        if index_path.is_file():
            with (output / "README.md").open("a", encoding="utf-8") as stream:
                for item in json.loads(index_path.read_text())["plots"]:
                    if item["status"] == "generated":
                        stream.write(f"\n![{item['name']}]({item['name']}.png)\n")
                    else:
                        stream.write(f"\n{item['name']}: not plotted ({item['reason']}).\n")
    else:
        records.save(
            output / "plot_index.json",
            {
                "schema": "intrance_run_plots_v1",
                "plots": [
                    {
                        "name": name,
                        "status": "skipped",
                        "files": [],
                        "reason": "No retained trajectory or local result",
                    }
                    for name in (
                        "trajectory_quality",
                        "controls",
                        "refinement_progress",
                        "numerical_checks",
                    )
                ],
                "sources": {},
                "new_propagation": False,
                "new_solver_or_native_run": False,
            },
        )


def _run_file(work, relative):
    """Use files inside this run, never an old absolute path stored in a receipt."""
    relative = Path(relative)
    if relative.is_absolute() or ".." in relative.parts:
        raise ValueError("Plot input must be local to the run")
    path = work
    for part in relative.parts:
        path = path / part
        if path.is_symlink():
            raise ValueError("Linked plot input: " + str(relative))
    if not path.is_file():
        raise FileNotFoundError(str(relative))
    return path


def _run_json(work, relative, sources):
    path = _run_file(work, relative)
    raw = path.read_bytes()
    sources[str(relative)] = {"sha256": hashlib.sha256(raw).hexdigest(), "bytes": len(raw)}
    value = json.loads(raw, object_pairs_hook=records._pairs)
    if not isinstance(value, dict):
        raise ValueError("Expected JSON object: " + str(relative))
    return value


def _run_table(work, relative, fields, sources):
    path = _run_file(work, relative)
    current = {}
    value = read_table(path, fields, current)
    sources[str(relative)] = current[path.name]
    return value


def _stage_key(value):
    if not isinstance(value, str) or not re.fullmatch(r"restore|quality_\d{2,3}(?:_\d{2})?", value):
        raise ValueError("Invalid selected stage key")
    return value


def _native_plot_history(work, label, sources):
    if label not in ("baseline", "candidate"):
        raise ValueError("Unknown native history")
    directory = f"native/stages/precise_{label}/result/"
    fields = ("MJD_d", "x_AU", "y_AU", "z_AU", "vx_AU_per_d", "vy_AU_per_d", "vz_AU_per_d")
    states = _run_table(work, directory + "bridge_states.csv", fields, sources)
    if len(states) < 2 or np.any(np.diff(states[:, 0]) <= 0):
        raise ValueError("Native state epochs must strictly increase")
    return {
        "states": states[:, 1:],
        "epochs": states[:, 0],
        "time": states[:, 0] - states[0, 0],
        "directory": directory,
    }


def _selected_plot_candidate(work, selection, sources):
    if not selection.get("selected_stage_key"):
        raise FileNotFoundError("No selected local stage was retained")
    key = _stage_key(selection["selected_stage_key"])
    result = _run_json(work, f"local_{key}/result.json", sources)
    candidate = result.get("candidate")
    if not isinstance(candidate, dict):
        raise ValueError("Selected stage has no returned candidate")
    return key, candidate


def _check_candidate_time(candidate, history):
    duration = float(candidate["duration_days"])
    bounds = np.asarray(candidate["duration_bounds_days"], dtype=float)
    if (
        bounds.shape != (2,)
        or not np.isfinite(bounds).all()
        or bounds[0] > bounds[1]
        or not bounds[0] <= duration <= bounds[1]
    ):
        raise ValueError("Candidate duration is outside its recorded bounds")
    epochs = history["epochs"]
    expected = candidate["departure_MJD"] + duration * np.arange(len(epochs)) / (len(epochs) - 1)
    if not np.allclose(epochs, expected, rtol=0, atol=1e-9):
        raise ValueError("Candidate and native history timing differ")


def _control_plot_data(work, selection, sources):
    history = _native_plot_history(work, "candidate", sources)
    commands = _run_table(
        work,
        history["directory"] + "bridge_commands.csv",
        ("MJD_d", "clock_rad", "cone_rad"),
        sources,
    )
    if commands.shape[0] != len(history["epochs"]) or not np.allclose(
        commands[:, 0], history["epochs"], rtol=0, atol=1e-9
    ):
        raise ValueError("Native control and state timing differ")
    angles = commands[:-1, 1:]
    if np.any(np.abs(angles[:, 0]) > math.pi + 1e-12) or np.any(
        (angles[:, 1] < -1e-12) | (angles[:, 1] > math.pi / 2 + 1e-12)
    ):
        raise ValueError("Native sail angles lie outside their bounds")
    _, candidate = _selected_plot_candidate(work, selection, sources)
    _check_candidate_time(candidate, history)
    saved = np.asarray(candidate["controls_clock_cone_rad"], dtype=float)
    if saved.shape != angles.shape or not np.allclose(saved, angles, rtol=0, atol=1e-12):
        raise ValueError("Selected controls differ from native commands")
    alpha, beta = angles.T
    normals = np.column_stack(
        (np.cos(beta), np.sin(beta) * np.cos(alpha), np.sin(beta) * np.sin(alpha))
    )
    turns = np.rad2deg(np.arccos(np.clip(np.sum(normals[1:] * normals[:-1], axis=1), -1, 1)))
    policy = candidate.get("control_policy")
    limit = None
    if policy and policy.get("mode") == "reference_relative_command_grid":
        duration = float(history["time"][-1])
        limit = math.degrees(
            min(
                float(policy["max_step_rad"]),
                float(policy["max_step_rate_rad_per_day"]) * duration / len(angles),
            )
        )
        if not math.isfinite(limit) or limit < 0:
            raise ValueError("Invalid saved normal-step bound")
    return {
        **history,
        "angles_deg": np.rad2deg(angles),
        "normal_steps_deg": turns,
        "allowed_step_deg": limit,
        "holds": len(angles),
    }


def _iteration_trace(text, phase, time_scale=None):
    """These are printed iterates, not all line-search trial points."""
    rows = []
    for line in text.splitlines():
        match = re.match(r"^\s*(\d+)(r?)\s+([-+\d.eE]+)\s+([-+\d.eE]+)\s+([-+\d.eE]+)\s", line)
        if not match:
            continue
        values = [float(match[i]) for i in (3, 4, 5)]
        if not all(math.isfinite(v) for v in values) or min(values[1:]) < 0:
            raise ValueError("Nonfinite or negative iteration residual")
        row = dict(
            iteration=int(match[1]),
            ipopt_restoration=bool(match[2]),
            objective=values[0],
            inf_pr=values[1],
            inf_du=values[2],
        )
        if phase == "quality":
            if time_scale is None or not math.isfinite(time_scale) or time_scale <= 0:
                raise ValueError("Time-stage trace lacks its objective scale")
            row["approx_duration_days"] = values[0] * time_scale
        rows.append(row)
    return rows


def _refinement_plot_data(work, selection, sources):
    points, notes = [], []
    keys = ["restore"]
    for item in selection.get("cells", []):
        if "stage_key" in item:
            key = _stage_key(item["stage_key"])
            if key not in keys:
                keys.append(key)
        else:
            continue
        if "selected_duration_days" not in item:
            continue
        selected = _stage_key(item.get("selected_stage_key", key))
        try:
            bundle = _run_json(work, f"local_{selected}/result.json", sources)
            candidate = bundle["candidate"]
            duration = float(item["selected_duration_days"])
            attempt_input = _run_json(work, f"local_{key}/input_candidate.json", sources)
            bounds = np.asarray(attempt_input["duration_bounds_days"], dtype=float)
            if (
                not math.isfinite(duration)
                or bounds.shape != (2,)
                or not np.isfinite(bounds).all()
                or bounds[0] > bounds[1]
                or not bounds[0] - 1e-7 <= float(attempt_input["duration_days"]) <= bounds[1] + 1e-7
                or abs(duration - float(candidate["duration_days"])) > 1e-7
            ):
                raise ValueError("Retained duration/bounds differ from the saved attempt")
            points.append(
                {
                    "stage": key,
                    "selected_stage": selected,
                    "duration_days": duration,
                    "bounds_days": bounds.tolist(),
                }
            )
        except FileNotFoundError as error:
            notes.append(f"Retained point unavailable: {error}")
    if selection.get("selected_stage_key"):
        try:
            restore = _run_json(work, "local_restore/result.json", sources).get("candidate")
            if restore is not None:
                duration = float(restore["duration_days"])
                if not math.isfinite(duration):
                    raise ValueError("Nonfinite restoration duration")
                points.insert(
                    0,
                    {
                        "stage": "restore",
                        "selected_stage": "restore",
                        "duration_days": duration,
                        "bounds_days": [duration, duration],
                    },
                )
        except FileNotFoundError:
            notes.append("Restoration point was not retained")
    # An interrupted worker can leave a real log before a selection ledger is written.
    attempted = []
    for directory in work.glob("local_quality_*"):
        if directory.is_dir() and (directory / "stdout.txt").is_file():
            key = _stage_key(directory.name.removeprefix("local_"))
            parts = key.split("_")
            attempted.append(((int(parts[1]), int(parts[2]) if len(parts) > 2 else 0), key))
    for _, key in sorted(attempted):
        if key not in keys:
            keys.append(key)
    keys = ["restore"] + sorted(
        set(keys) - {"restore"},
        key=lambda key: (
            int(key.split("_")[1]),
            int(key.split("_")[2]) if len(key.split("_")) > 2 else 0,
        ),
    )
    trace, trace_key, phase = [], None, None
    for key in reversed(keys):
        try:
            log = _run_file(work, f"local_{key}/stdout.txt")
            raw = log.read_bytes()
            sources[f"local_{key}/stdout.txt"] = {
                "sha256": hashlib.sha256(raw).hexdigest(),
                "bytes": len(raw),
            }
            try:
                bundle = _run_json(work, f"local_{key}/result.json", sources)
                phase = bundle.get("result", {}).get("phase")
                scale = bundle.get("attempt", {}).get("protocol", {}).get("time_scale_days")
            except FileNotFoundError:
                context = _run_json(work, f"local_{key}/context.json", sources)
                candidate = _run_json(work, f"local_{key}/input_candidate.json", sources)
                phase = context.get("phase")
                scale = objective.dynamical_time_days(candidate["mission_profile"])
            if phase not in ("restore", "quality"):
                raise ValueError("Saved solver phase is missing or unknown")
            trace = _iteration_trace(raw.decode("utf-8", errors="replace"), phase, scale)
            if trace:
                trace_key = key
                break
            notes.append("No scalar iteration rows in " + key)
        except FileNotFoundError:
            continue
    if not trace:
        notes.append("No retained scalar Ipopt trace")
    if not points and not trace:
        raise FileNotFoundError("No retained refinement points or scalar solver trace")
    return {
        "points": points,
        "trace": trace,
        "trace_stage": trace_key,
        "trace_phase": phase,
        "notes": notes,
    }


def _error_history(first, second):
    a, b = np.asarray(first, dtype=float), np.asarray(second, dtype=float)
    if (
        a.ndim != 2
        or b.ndim != 2
        or a.shape[0] != b.shape[0]
        or min(a.shape[1], b.shape[1]) < 6
        or not np.isfinite(a).all()
        or not np.isfinite(b).all()
    ):
        raise ValueError("Numerical histories have different shapes or nonfinite values")
    difference = a[:, :6] - b[:, :6]
    return {
        "position_m": np.linalg.norm(difference[:, :3], axis=1) * AU_M,
        "velocity_m_s": np.linalg.norm(difference[:, 3:6], axis=1) * AU_M / DAY_S,
    }


def _numerical_plot_data(work, selection, sources):
    review = _run_json(work, "candidate_review/native_comparison.json", sources)
    contract = review["validation_contract"]
    definitions = [
        (
            "native",
            "Native − RK4-256",
            review.get("native_vs_continuous", {}).get("256"),
            [contract["native_agreement_position_m"], contract["native_agreement_velocity_m_s"]],
        ),
        (
            "mesh",
            "RK4-128 − RK4-256",
            review.get("fine_mesh128_vs256"),
            [contract["fine_mesh_position_m"], contract["fine_mesh_velocity_m_s"]],
        ),
        (
            "nodes",
            "Nodes − RK4-256",
            review.get("nodes_vs_continuous256"),
            [100.0, 0.001] if review.get("node_continuity_gate_applied") is True else [None, None],
        ),
    ]
    histories, missing = {}, []
    for mesh in (128, 256):
        relative = f"candidate_review/continuous_RK{mesh}.npy"
        try:
            raw = _run_file(work, relative).read_bytes()
            digest = hashlib.sha256(raw).hexdigest()
            binding = review.get("saved_histories", {}).get(str(mesh), {})
            if binding.get("sha256") != digest:
                raise ValueError("Saved history hash differs: " + relative)
            values = np.load(io.BytesIO(raw), allow_pickle=False)
            if values.ndim != 2 or values.shape[1] < 6 or not np.isfinite(values).all():
                raise ValueError("Invalid saved continuous history")
            histories[mesh] = values
            sources[relative] = {"sha256": digest, "bytes": len(raw)}
        except FileNotFoundError:
            missing.append(relative)
    native, candidate = None, None
    try:
        native = _native_plot_history(work, "candidate", sources)
        _, candidate = _selected_plot_candidate(work, selection, sources)
        _check_candidate_time(candidate, native)
    except FileNotFoundError as error:
        missing.append(str(error))
    if native is not None:
        time_days = native["time"]
        for values in histories.values():
            if len(values) != len(time_days):
                raise ValueError("Continuous and native history counts differ")
    else:
        time_days = None
    comparisons = []
    for key, label, recorded, gates in definitions:
        series = None
        if time_days is not None and 256 in histories:
            if key == "native":
                series = _error_history(native["states"], histories[256])
            elif key == "mesh" and 128 in histories:
                series = _error_history(histories[128], histories[256])
            elif key == "nodes" and candidate is not None:
                series = _error_history(candidate["nodes_cartesian"], histories[256])
        maxima = (
            {name: float(np.max(values)) for name, values in series.items()} if series else recorded
        )
        if maxima is None:
            comparisons.append({"key": key, "label": label, "mode": "unavailable", "gates": gates})
            continue
        for name in ("position_m", "velocity_m_s"):
            if not math.isfinite(maxima[name]) or maxima[name] < 0:
                raise ValueError("Invalid recorded numerical maximum")
            if (
                series is not None
                and recorded is not None
                and not math.isclose(
                    maxima[name],
                    recorded[name],
                    rel_tol=1e-8,
                    abs_tol=1e-6 if name == "position_m" else 1e-11,
                )
            ):
                raise ValueError("Saved numerical maximum differs from its history: " + key)
        comparisons.append(
            {
                "key": key,
                "label": label,
                "mode": "history" if series is not None else "recorded_maximum",
                "series": series,
                "maxima": maxima,
                "gates": gates,
            }
        )
    if all(row["mode"] == "unavailable" for row in comparisons):
        raise FileNotFoundError("No numerical histories or recorded maxima")
    return {
        "time_days": time_days,
        "comparisons": comparisons,
        "missing": missing,
        "scaled_gap_limit": contract.get("continuity_scaled_limit"),
    }


def plot(workspace, output=None):
    """Render retained data only. An output directory leaves an old run unchanged."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    work = Path(workspace).resolve(strict=True)
    explicit_output = output is not None
    output = Path(output) if explicit_output else work / "report"
    output.mkdir(parents=True, exist_ok=not explicit_output)
    sources, plots = {}, []
    try:
        data = _run_json(work, "report/metrics_before_timing.json", sources)
    except FileNotFoundError:
        data = {"trajectories": []}
    selection = data.get("selection") or {}
    if not selection and (work / "selection.json").is_file():
        selection = _run_json(work, "selection.json", sources)
    plt.rcParams.update(
        {
            "font.size": 10,
            "axes.titlesize": 11,
            "axes.labelsize": 10,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "pdf.fonttype": 42,
            "svg.fonttype": "none",
        }
    )

    def save(fig, name, definition, details=None):
        files = []
        for extension in ("png", "svg", "pdf"):
            path = output / f"{name}.{extension}"
            fig.savefig(path, dpi=180, bbox_inches="tight")
            files.append(path.name)
        plt.close(fig)
        plots.append(
            {
                "name": name,
                "status": "generated",
                "files": files,
                "definition": definition,
                **(details or {}),
            }
        )

    def skip(name, error):
        plots.append({"name": name, "status": "skipped", "files": [], "reason": str(error)})

    try:
        rows = [row for row in data.get("trajectories", []) if row["label"] != "native_selected"]
        if not rows:
            raise FileNotFoundError("No retained native trajectory rows")
        fig, axes = plt.subplots(
            1, 3, figsize=(10.5, 3.9), gridspec_kw={"width_ratios": [1.3, 1, 1]}
        )
        palette = {"baseline": "#b35b30", "candidate": "#168c8c"}
        for row in rows:
            history = _native_plot_history(work, row["label"], sources)
            if abs(history["time"][-1] - row["metrics"]["flight_time_days"]) > 1e-7:
                raise ValueError("Native history and reported duration differ")
            x, y = history["states"][:, :2].T
            axes[0].plot(x, y, label=row["label"], color=palette[row["label"]], lw=1.3)
            axes[0].scatter(x[-1], y[-1], color=palette[row["label"]], s=15)
        axes[0].scatter(0, 0, marker="*", color="#d7a316", s=60, label="Sun")
        axes[0].set(xlabel="x [AU]", ylabel="y [AU]", title="Native path", aspect="equal")
        axes[0].legend(fontsize=8)
        for axis, key, title in [
            (axes[1], "native_fitness", "Native fitness"),
            (axes[2], "own_fitness", "Time score"),
        ]:
            for i, row in enumerate(rows):
                value = row.get(key)
                if value is None:
                    axis.text(i, 0.02, "not available", rotation=90, ha="center")
                else:
                    axis.bar(i, value, color=palette[row["label"]], width=0.55)
                    axis.text(i, value, f"{value:.4g}", ha="center", va="bottom", fontsize=9)
            axis.set_xticks(range(len(rows)), [row["label"] for row in rows])
            axis.set_title(title)
            axis.margins(y=0.2)
        fig.tight_layout()
        save(fig, "trajectory_quality", "Native paths and two distinct reporting scores.")
    except (FileNotFoundError, ValueError, KeyError) as error:
        plt.close("all")
        skip("trajectory_quality", error)

    try:
        controls = _control_plot_data(work, selection, sources)
        fig, axes = plt.subplots(3, 1, figsize=(7.2, 5.5), sharex=True)
        for i, label in enumerate(("Clock [deg]", "Cone [deg]")):
            values = controls["angles_deg"][:, i]
            axes[i].step(controls["time"], np.r_[values, values[-1]], where="post", lw=1.1)
            axes[i].set_ylabel(label)
        axes[2].plot(controls["time"][1:-1], controls["normal_steps_deg"], lw=1.1)
        if controls["allowed_step_deg"] is not None:
            axes[2].axhline(
                controls["allowed_step_deg"],
                ls="--",
                color="#777777",
                label=f"Bound {controls['allowed_step_deg']:.4f} deg",
            )
            axes[2].legend(fontsize=8)
        axes[2].set(xlabel="Flight time [days]", ylabel="Normal step [deg]")
        axes[0].set_title(f"Commands: {controls['holds']} holds")
        for axis in axes:
            axis.grid(alpha=0.2)
        fig.tight_layout()
        save(
            fig,
            "controls",
            "Held RTN commands; normal steps are not hardware slew rates.",
            {"holds": controls["holds"], "allowed_normal_step_deg": controls["allowed_step_deg"]},
        )
    except (FileNotFoundError, ValueError, KeyError) as error:
        plt.close("all")
        skip("controls", error)

    try:
        refinement = _refinement_plot_data(work, selection, sources)
        fig = plt.figure(figsize=(7.2, 5.5))
        grid = fig.add_gridspec(2, 2)
        stages, value_axis, residual = (
            fig.add_subplot(grid[0, :]),
            fig.add_subplot(grid[1, 0]),
            fig.add_subplot(grid[1, 1]),
        )
        points, trace = refinement["points"], refinement["trace"]
        if points:
            for i, point in enumerate(points):
                stages.vlines(i, *point["bounds_days"], color="#bbbbbb", linewidth=10, alpha=0.5)
            stages.plot(range(len(points)), [p["duration_days"] for p in points], "o-", ms=4)
            stages.set_xticks(
                range(len(points)),
                ["Restore" if p["stage"] == "restore" else str(i) for i, p in enumerate(points)],
            )
            stages.annotate(
                f"{points[-1]['duration_days']:.3f}",
                (len(points) - 1, points[-1]["duration_days"]),
                xytext=(0, 7),
                textcoords="offset points",
                ha="center",
                fontsize=9,
            )
        else:
            stages.text(0.5, 0.5, "No retained selection", ha="center", transform=stages.transAxes)
        stages.set(
            title="Retained time and cell bounds", ylabel="TOF [days]", xlabel="Solve attempt"
        )
        if trace:
            iterations = [r["iteration"] for r in trace]
            time_phase = refinement["trace_phase"] == "quality"
            field = "approx_duration_days" if time_phase else "objective"
            value_axis.plot(iterations, [r[field] for r in trace])
            value_axis.set(ylabel="Approx. TOF [days]" if time_phase else "Arrival objective")
            inf_pr = [r["inf_pr"] for r in trace]
            residual.plot(iterations, inf_pr)
            if max(inf_pr) > 0:
                residual.set_yscale("symlog", linthresh=1e-15)
                residual.yaxis.get_major_locator().set_params(numticks=7)
            residual.set_ylabel("NLP inf_pr")
            parts = refinement["trace_stage"].split("_")
            label = (
                "Restoration"
                if parts[0] == "restore"
                else f"Time solve {int(parts[1])}.{int(parts[2]) if len(parts) > 2 else 1}"
            )
            value_axis.set_title(label, fontsize=9)
        else:
            for axis in (value_axis, residual):
                axis.text(0.5, 0.5, "No saved iteration log", ha="center", transform=axis.transAxes)
        for axis in (value_axis, residual):
            axis.set_xlabel("Ipopt iteration")
        for axis in (stages, value_axis, residual):
            axis.grid(alpha=0.2)
        fig.tight_layout()
        save(
            fig,
            "refinement_progress",
            "Retained points and the last saved scalar Ipopt trace.",
            {
                "selected_stage": selection.get("selected_stage_key"),
                "trace_stage": refinement["trace_stage"],
                "trace_phase": refinement["trace_phase"],
                "notes": refinement["notes"],
                "iterate_tof_approximate_from_printed_objective": refinement["trace_phase"]
                == "quality",
            },
        )
    except (FileNotFoundError, ValueError, KeyError) as error:
        plt.close("all")
        skip("refinement_progress", error)

    try:
        numerical = _numerical_plot_data(work, selection, sources)
        fig, axes = plt.subplots(3, 2, figsize=(7.2, 7.0))
        summaries = []
        for row, comparison in enumerate(numerical["comparisons"]):
            for column, (field, unit) in enumerate((("position_m", "m"), ("velocity_m_s", "m/s"))):
                axis = axes[row, column]
                axis.set_title(comparison["label"])
                axis.set_ylabel(f"Position [{unit}]" if column == 0 else f"Velocity [{unit}]")
                if comparison["mode"] == "unavailable":
                    axis.text(0.5, 0.5, "Not retained", ha="center", transform=axis.transAxes)
                    continue
                gate = comparison["gates"][column]
                if comparison["mode"] == "history":
                    values = comparison["series"][field]
                    axis.plot(numerical["time_days"], values, lw=1.0)
                    axis.set_xlabel("Flight time [days]")
                else:
                    axis.bar([0], [comparison["maxima"][field]], width=0.4)
                    axis.set_xticks([0], ["Saved max"])
                if gate is not None:
                    axis.axhline(gate, ls="--", color="#777777", label=f"Limit {gate:g}")
                    axis.legend(fontsize=8)
                if comparison["maxima"][field] > 0:
                    axis.set_yscale("symlog", linthresh=1e-6 if column == 0 else 1e-12)
                axis.grid(alpha=0.2)
            summaries.append({key: value for key, value in comparison.items() if key != "series"})
        fig.tight_layout()
        save(
            fig,
            "numerical_checks",
            "Three separate physical-unit numerical comparisons.",
            {
                "comparisons": summaries,
                "missing_sources": numerical["missing"],
                "scaled_continuity_limit_separate_from_physical_curves": numerical[
                    "scaled_gap_limit"
                ],
            },
        )
    except (FileNotFoundError, ValueError, KeyError) as error:
        plt.close("all")
        skip("numerical_checks", error)
    index = {
        "schema": "intrance_run_plots_v1",
        "plots": plots,
        "sources": sources,
        "new_propagation": False,
        "new_solver_or_native_run": False,
    }
    (output / "plot_index.json").write_text(json.dumps(index, indent=2, allow_nan=False) + "\n")
    return index


import hashlib
import io
import numpy as np

AU_M = 149597870691.0
DAY_S = 86400.0


def read_table(path, columns, sources):
    raw = path.read_bytes()
    sources[path.name] = {"sha256": hashlib.sha256(raw).hexdigest(), "bytes": len(raw)}
    reader = csv.DictReader(io.StringIO(raw.decode("utf-8-sig")), skipinitialspace=True)
    names = reader.fieldnames or []
    if len(names) != len(set(names)) or not set(columns).issubset(names):
        raise ValueError(f"Missing or duplicate CSV columns: {path}")
    rows = []
    for row in reader:
        if None in row or any((row.get(name) is None for name in names)):
            raise ValueError(f"Incomplete CSV row: {path}")
        rows.append([float(row[name]) for name in columns])
    values = np.asarray(rows, dtype=float)
    if values.ndim != 2 or not len(values) or (not np.isfinite(values).all()):
        raise ValueError(f"Empty or nonfinite CSV: {path}")
    return values


def load_saved_result(directory, metrics_path=None):
    """Read and bind the saved state, command and final-target records once."""
    directory = Path(directory)
    sources = {}
    state_columns = ("MJD_d", "x_AU", "y_AU", "z_AU", "vx_AU_per_d", "vy_AU_per_d", "vz_AU_per_d")
    states = read_table(directory / "native/bridge_states.csv", state_columns, sources)
    commands = read_table(
        directory / "native/bridge_commands.csv", ("MJD_d", "clock_rad", "cone_rad"), sources
    )
    endpoint_columns = (
        "requested_final_MJD_d",
        "actual_final_MJD_d",
        "target_MJD_d",
        "target_body_id",
        "spacecraft_x_AU",
        "spacecraft_y_AU",
        "spacecraft_z_AU",
        "spacecraft_vx_AU_per_d",
        "spacecraft_vy_AU_per_d",
        "spacecraft_vz_AU_per_d",
        "target_x_AU",
        "target_y_AU",
        "target_z_AU",
        "target_vx_AU_per_d",
        "target_vy_AU_per_d",
        "target_vz_AU_per_d",
    )
    endpoints = read_table(
        directory / "native/bridge_target_endpoint.csv", endpoint_columns, sources
    )
    if len(states) < 2 or len(states) != len(commands) or len(endpoints) != 1:
        raise ValueError("Expected aligned state/command boundaries and one target endpoint")
    epochs = states[:, 0]
    if np.any(np.diff(epochs) <= 0) or not np.allclose(epochs, commands[:, 0], rtol=0, atol=1e-09):
        raise ValueError("State and command epochs must align and strictly increase")
    endpoint = endpoints[0]
    if endpoint[3] != 2 or not np.allclose(endpoint[:3], epochs[-1], rtol=0, atol=1e-09):
        raise ValueError("Mercury endpoint epochs differ from the final boundary")
    if not np.allclose(states[-1, 1:], endpoint[4:10], rtol=0, atol=1e-12):
        raise ValueError("Endpoint spacecraft state differs from the recorded history")
    time = epochs - epochs[0]
    duration = float(time[-1])
    angles = commands[:-1, 1:]
    clock, cone = (angles[:, 0], angles[:, 1])
    if np.any(cone < -1e-12) or np.any(cone > np.pi / 2 + 1e-12):
        raise ValueError("Cone commands lie outside the ideal-sail 0..pi/2 range")
    normals = np.column_stack(
        (np.cos(cone), np.sin(cone) * np.cos(clock), np.sin(cone) * np.sin(clock))
    )
    turns = np.rad2deg(np.arccos(np.clip(np.sum(normals[1:] * normals[:-1], axis=1), -1, 1)))
    position_km = float(np.linalg.norm(endpoint[4:7] - endpoint[10:13]) * AU_M / 1000)
    velocity_m_s = float(np.linalg.norm(endpoint[7:10] - endpoint[13:16]) * AU_M / DAY_S)
    summary = {
        "schema": "saved_schedule_plot_v1",
        "flight_time_days": duration,
        "departure_MJD": float(epochs[0]),
        "boundary_records": len(states),
        "held_commands": len(angles),
        "arrival_position_km": position_km,
        "arrival_velocity_m_s": velocity_m_s,
        "maximum_same_frame_normal_step_deg": float(turns.max()) if len(turns) else None,
        "units": {
            "position": "AU",
            "time": "day",
            "input_angles": "radian",
            "plot_angles": "degree",
        },
        "new_propagation_or_optimisation": False,
        "new_native_validation": False,
        "command_steps_are_not_physical_slew": True,
        "source_files": sources,
    }
    metrics_path = Path(metrics_path) if metrics_path is not None else directory / "metrics.json"
    guard_limit = None
    if metrics_path.exists():
        raw = metrics_path.read_bytes()
        metrics = json.loads(raw)
        candidate = metrics["candidate"]
        measured = (duration, position_km, velocity_m_s)
        declared = [
            candidate[key] for key in ("duration_days", "position_miss_km", "velocity_miss_m_s")
        ]
        if not np.isfinite(declared).all() or not np.allclose(
            measured, declared, rtol=0, atol=(1e-07, 1e-06, 1e-08)
        ):
            raise ValueError("Metrics disagree with the saved endpoint or flight time")
        sources["metrics.json"] = {"sha256": hashlib.sha256(raw).hexdigest(), "bytes": len(raw)}
        guard = metrics.get("control_guard", {})
        if "allowed_normal_step_deg" in guard:
            guard_limit = float(guard["allowed_normal_step_deg"])
            if not np.isfinite(guard_limit) or guard_limit < 0:
                raise ValueError("Invalid recorded command-step limit")
    return {
        "states": states,
        "time": time,
        "angles": angles,
        "normal_steps": turns,
        "target_r_AU": endpoint[10:13],
        "guard_limit_deg": guard_limit,
        "summary": summary,
    }


def _saved_trajectory_figure(data):
    """Draw the stored XY path and its departure/arrival markers."""
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MultipleLocator

    position = data["states"][:, 1:4]
    target = data["target_r_AU"]
    duration = data["summary"]["flight_time_days"]
    blue, orange = "#0072B2", "#D55E00"
    fig, ax = plt.subplots(figsize=(7.2, 7.2))
    ax.plot(position[:, 0], position[:, 1], color=blue, lw=1.4, label="Spacecraft")
    ax.scatter(
        *position[0, :2], color="#222222", marker="s", s=34, label="Earth departure", zorder=4
    )
    ax.scatter(*target[:2], color=orange, marker="x", s=50, label="Mercury arrival", zorder=5)
    ax.scatter(
        *position[-1, :2],
        edgecolors=blue,
        facecolors="none",
        s=72,
        label="Spacecraft arrival",
        zorder=6,
    )
    ax.scatter(0, 0, color="#D89C12", marker="*", s=95, label="Sun", zorder=4)
    extent = float(np.max(np.abs(np.r_[position[:, :2].ravel(), target[:2]])))
    limit = max(0.25, math.ceil(1.08 * extent * 4) / 4)
    ax.set(
        xlim=(-limit, limit),
        ylim=(-limit, limit),
        xlabel="x [AU]",
        ylabel="y [AU]",
        title=f"Earth-Mercury\nTOF {duration:.3f} days",
    )
    ax.set_aspect("equal", adjustable="box")
    ax.xaxis.set_major_locator(MultipleLocator(0.5))
    ax.yaxis.set_major_locator(MultipleLocator(0.5))
    ax.grid(color="#d5d5d5", alpha=0.7, lw=0.6)
    ax.set_facecolor("white")
    ax.legend(loc="upper right", fontsize=8.5, frameon=False)
    ax.tick_params(labelsize=10)
    ax.xaxis.label.set_size(11)
    ax.yaxis.label.set_size(11)
    ax.title.set_size(13)
    fig.subplots_adjust(left=0.13, right=0.96, bottom=0.09, top=0.92)
    return fig


def draw_saved_result(data, output):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    output = Path(output)
    output.mkdir(parents=True, exist_ok=False)
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "savefig.facecolor": "white",
        }
    )
    summary = data["summary"]
    duration = summary["flight_time_days"]
    files = {}

    def save(fig, name):
        for suffix in ("png", "svg", "pdf"):
            path = output / f"{name}.{suffix}"
            fig.savefig(path, dpi=300, bbox_inches="tight")
            files[path.name] = hashlib.sha256(path.read_bytes()).hexdigest()
        plt.close(fig)

    blue, orange = ("#0072B2", "#D55E00")
    save(_saved_trajectory_figure(data), "trajectory")
    fig, axes = plt.subplots(3, 1, figsize=(8, 7), sharex=True)
    angles = np.rad2deg(data["angles"])
    for index, label in enumerate(("Clock [deg]", "Cone [deg]")):
        axes[index].step(
            data["time"],
            np.r_[angles[:, index], angles[-1, index]],
            where="post",
            color=blue if index == 0 else orange,
            lw=1.1,
        )
        axes[index].set_ylabel(label)
    axes[2].plot(
        data["time"][1:-1],
        data["normal_steps"],
        color="#009E73",
        lw=1.1,
        label="Adjacent command normals",
    )
    if data["guard_limit_deg"] is not None:
        axes[2].axhline(
            data["guard_limit_deg"],
            color="#555555",
            ls="--",
            lw=0.8,
            label="Recorded command-grid bound",
        )
    axes[2].set(xlabel="Elapsed flight time [days]", ylabel="Normal step [deg]", xlim=(0, duration))
    axes[2].legend(loc="upper right", fontsize=8)
    for ax in axes:
        ax.grid(color="#cccccc", alpha=0.5, lw=0.6)
    fig.suptitle(f"Commands: {summary['held_commands']} holds", fontsize=13)
    fig.tight_layout(rect=(0, 0.065, 1, 0.97))
    save(fig, "controls")
    summary.update(
        output_sha256=files,
        plot_environment={"numpy": np.__version__, "matplotlib": matplotlib.__version__},
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    )
    internal = output / "_internal"
    internal.mkdir(exist_ok=True)
    (internal / "plot_checks.json").write_text(
        json.dumps(summary, indent=2, allow_nan=False) + "\n"
    )
    public = {
        key: value
        for key, value in summary.items()
        if key not in {"source_files", "output_sha256", "script_sha256"}
    }
    public["output_files"] = list(files)
    (output / "plot_summary.json").write_text(json.dumps(public, indent=2, allow_nan=False) + "\n")
    return public
