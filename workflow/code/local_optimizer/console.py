"""Readable progress for the console; these messages never decide acceptance."""

import codecs
import json
import math
from pathlib import Path
import re
import subprocess
import sys
import time


class TextTail:
    """Read complete new lines while a child is still writing its log."""

    def __init__(self, path):
        self.path = Path(path)
        self.offset = 0
        self.pending = ""
        self.decoder = codecs.getincrementaldecoder("utf-8")(errors="replace")

    def read(self, *, final=False):
        try:
            with self.path.open("rb") as stream:
                if stream.seek(0, 2) < self.offset:
                    self.offset, self.pending = 0, ""
                    self.decoder.reset()
                stream.seek(self.offset)
                data = stream.read(65536)
                self.offset = stream.tell()
                at_end = self.offset == stream.seek(0, 2)
        except OSError:
            # A temporarily unavailable progress log must not stop the calculation.
            return []
        text = self.pending + self.decoder.decode(data)
        parts = text.split("\n")
        self.pending = parts.pop()
        if final and at_end and self.pending:
            parts.append(self.pending)
            self.pending = ""
        return [line.rstrip("\r") for line in parts]


class Console:
    def __init__(self, workspace):
        self.workspace = Path(workspace)
        self.started = time.monotonic()
        self.log_warning_shown = False

    def emit(self, message, *, mirror=True, preserve_spacing=False):
        # Plain text is also read by the shared Windows window, never as a command.
        text = str(message) if preserve_spacing else " ".join(str(message).split())
        message = re.sub(r"[^\x20-\x7e]", "?", text)
        elapsed = int(time.monotonic() - self.started)
        line = f"[{elapsed // 60:02d}:{elapsed % 60:02d}] {message}"
        print(line, flush=True)
        if mirror:
            try:
                with (self.workspace / "progress.log").open("a", encoding="ascii") as stream:
                    stream.write(line + "\n")
            except OSError:
                if not self.log_warning_shown:
                    print(
                        (
                            "Progress could not be copied to the other window; this "
                            "console still shows the calculation."
                        ),
                        file=sys.stderr,
                        flush=True,
                    )
                    self.log_warning_shown = True

    def setup(self, config):
        self.emit("Earth -> Mercury | ideal solar sail")
        self.emit(
            f"ENC seed {config['seed']} | stop after {config['evaluation_limit']} evaluations | "
            f"export {config['horizon_days']:g} days and {config['holds']} held commands"
        )
        mission = config["mission_profile"]
        self.emit(
            f"Local arrival limits: {mission['position_target_m'] / 1000:,.6f} km | "
            f"{mission['velocity_target_m_s']:.9f} m/s"
        )

    def local_result(self, name, bundle):
        result, candidate, review = (
            bundle["result"],
            bundle.get("candidate"),
            bundle.get("host_validation"),
        )
        self.emit(
            f"{name}: {result.get('solver_status', 'no solver status')} | "
            f"{result.get('solver_iterations', '?')} iterations"
        )
        if candidate is not None:
            self.emit(f"Returned TOF: {candidate['duration_days']:.9f} days")
        if review is not None:
            encounter = review.get("encounter", {})
            if "position_miss_km" in encounter and "velocity_miss_m_s" in encounter:
                self.emit(
                    f"Arrival from continuous host propagation: "
                    f"{encounter['position_miss_km']:,.6f} km | "
                    f"{encounter['velocity_miss_m_s']:.9f} m/s"
                )
            passed = review.get("accepted_for_native_review") is True
            self.emit(
                "Host checks passed; native replay still required."
                if passed
                else "Host checks did not pass; this candidate is not accepted."
            )


_ITERATION = re.compile(r"^\s*(\d+)(r?)\s+([-+\d.eE]+)\s+([-+\d.eE]+)\s+([-+\d.eE]+)\s+")


def solver_line(line, phase, time_scale_days):
    """Translate Ipopt's existing output, without adding a solver callback."""
    match = _ITERATION.match(line)
    if match is None:
        return None
    iteration, restoration, objective, primal, dual = match.groups()
    try:
        objective, primal, dual = map(float, (objective, primal, dual))
    except ValueError:
        return None
    if not all(math.isfinite(v) for v in (objective, primal, dual)):
        return None
    if phase == "quality":
        duration = objective * time_scale_days
        if not math.isfinite(duration):
            return None
        value = f"{duration:.6f}"
    else:
        value = f"{objective:.6g}"
    step = iteration + restoration
    return f"{step:>6} | {value:>15} | {primal:>15.2e} | {dual:>13.2e}"


def wait_for_solver(child, timeout, console, log_path, *, phase, time_scale_days, solver_limit):
    """Keep the existing process deadline while showing real child-log progress."""
    start = time.monotonic()
    tail = TextTail(log_path)
    latest, last_display, last_activity = None, start - 2, start
    header_printed = False

    def show_iteration(row):
        nonlocal header_printed
        if not header_printed:
            if phase == "quality":
                console.emit(
                    "TOF is the time being tried here. Full-flight checks follow this solve."
                )
                value_label = "TOF now [d]"
            else:
                console.emit("Fixed time first; arrival cost combines position and velocity miss.")
                value_label = "Arrival cost"
            console.emit("Constraint err.: mismatch in the model equations/limits.")
            console.emit("Dual err.: scaled stationarity check. r = Ipopt feasibility recovery.")
            console.emit(
                f"{'Step':>6} | {value_label:>15} | {'Constraint err.':>15} | {'Dual err.':>13}",
                preserve_spacing=True,
            )
            header_printed = True
        console.emit(row, preserve_spacing=True)

    while True:
        now = time.monotonic()
        for line in tail.read():
            value = solver_line(line, phase, time_scale_days)
            if value is not None:
                latest = value
            elif line.startswith("Local step:"):
                console.emit(line.removeprefix("Local step:").strip())
                last_activity = now
        if latest is not None and now - last_display >= 2:
            show_iteration(latest)
            latest, last_display, last_activity = None, now, now
        if now - last_activity >= 10:
            console.emit(
                f"Local stage running: {now - start:.0f} s elapsed | "
                f"solver limit {solver_limit:g} s; waiting for the next iteration output"
            )
            last_activity = now
        remaining = timeout - (time.monotonic() - start)
        if remaining <= 0:
            raise subprocess.TimeoutExpired(child.args, timeout)
        try:
            code = child.wait(timeout=min(0.25, remaining))
        except subprocess.TimeoutExpired:
            continue
        # Drain the last output before returning the actual process status.
        while True:
            lines = tail.read(final=True)
            if not lines:
                break
            for line in lines:
                value = solver_line(line, phase, time_scale_days)
                if value is not None:
                    latest = value
                elif line.startswith("Local step:"):
                    console.emit(line.removeprefix("Local step:").strip())
        if latest is not None:
            show_iteration(latest)
        return code


def _result_plots(workspace):
    report = Path(workspace) / "report"
    index_path = report / "plot_index.json"
    if report.is_symlink() or index_path.is_symlink():
        raise ValueError("Linked plot summary is not supported")
    names = set()
    pattern = re.compile(
        r"(trajectory_quality|controls|refinement_progress|numerical_checks)\.(png|svg|pdf)\Z"
    )
    if index_path.exists():
        index = json.loads(index_path.read_text(encoding="utf-8"))
        if (
            not isinstance(index, dict)
            or index.get("schema") != "intrance_run_plots_v1"
            or not isinstance(index.get("plots"), list)
        ):
            raise ValueError("Malformed plot index")
        for item in index["plots"]:
            if not isinstance(item, dict):
                raise ValueError("Malformed plot entry")
            if item.get("status") != "generated":
                continue
            files = item.get("files")
            if not isinstance(files, list):
                raise ValueError("Malformed generated plot list")
            for name in files:
                match = pattern.fullmatch(name) if isinstance(name, str) else None
                if match is None or match[1] != item.get("name"):
                    raise ValueError("Plot index contains an unexpected filename")
                names.add(name)
    else:
        # Older runs have no index, but their real plot files can still be listed.
        names = (
            {p.name for p in report.iterdir() if pattern.fullmatch(p.name)}
            if report.is_dir()
            else set()
        )
    return [
        str((report / name).relative_to(workspace))
        for name in sorted(names)
        if not (report / name).is_symlink() and (report / name).is_file()
    ]


def result_summary(completion, workspace):
    metrics = completion.get("six_metrics", {})
    mission = metrics.get("mission_profile", {})
    selected = metrics.get("selection") or {}
    plot_note = None
    try:
        plots = _result_plots(workspace)
    except (OSError, ValueError) as error:
        plots = []
        plot_note = "Plot list unavailable: " + str(error)
    result = {
        "status": "passed" if completion.get("workflow_complete") else "not completed",
        "scientific_checks_passed": completion.get("scientific_accepted") is True,
        "report_complete": completion.get("report_complete") is True,
        "flight_time_days": metrics.get("final_flight_time_days"),
        "arrival_position_miss_km": metrics.get("final_position_miss_km"),
        "arrival_velocity_miss_m_s": metrics.get("final_velocity_miss_m_s"),
        "arrival_position_limit_km": (
            mission["position_target_m"] / 1000 if "position_target_m" in mission else None
        ),
        "arrival_velocity_limit_m_s": mission.get("velocity_target_m_s"),
        "native_fitness": metrics.get("final_native_fitness"),
        "own_score": metrics.get("final_own_fitness"),
        "workflow_wall_seconds": completion.get("total_start_to_end_wall_seconds"),
        "local_solver_calls": selected.get("total_solver_calls"),
        "report": "report/README.md",
        "plots": plots,
    }
    if plot_note is not None:
        result["plot_list_note"] = plot_note
    return result


def finish(console, completion):
    summary = result_summary(completion, console.workspace)
    with (console.workspace / "result_summary.json").open("w", encoding="utf-8") as stream:
        json.dump(summary, stream, indent=2, allow_nan=False)
        stream.write("\n")
    console.emit(
        "Workflow passed."
        if completion.get("workflow_complete")
        else "Workflow did not complete successfully; read the saved report."
    )
    for key, label, unit in (
        ("flight_time_days", "Flight time", "days"),
        ("arrival_position_miss_km", "Arrival position miss", "km"),
        ("arrival_velocity_miss_m_s", "Arrival velocity miss", "m/s"),
        ("workflow_wall_seconds", "Complete workflow time (including waits)", "s"),
        ("native_fitness", "Native fitness", ""),
        ("own_score", "Own score", ""),
    ):
        value = summary[key]
        console.emit(
            f"{label}: {value:.9g} {unit}" if value is not None else f"{label}: unavailable"
        )
    if summary.get("plot_list_note"):
        console.emit(summary["plot_list_note"])
    console.emit("Outputs: result_summary.json and report/README.md")
    try:
        (console.workspace / "progress.done").write_text("done\n", encoding="ascii")
    except OSError:
        # The shared display may time out, but the saved numerical result is unchanged.
        print(
            "Final display notification unavailable; read result_summary.json for the status.",
            file=sys.stderr,
            flush=True,
        )


def print_entry(result):
    """Small CLI ending; the numerical result and raw records stay in files."""
    print(
        (
            "Run complete."
            if result.get("workflow_complete")
            else "Run stopped; inspect the saved report."
        ),
        flush=True,
    )
    if result.get("error"):
        print("Reason: " + str(result["error"]), flush=True)
