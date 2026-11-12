#!/usr/bin/env python3
"""Run or inspect the fixed Earth–Mercury shared calculation."""

import argparse
import json
from pathlib import Path
import sys

WORKFLOW_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(WORKFLOW_ROOT / "code"))


def main(argv=None):
    args = list(sys.argv[1:] if argv is None else argv)
    command = args.pop(0) if args and not args[0].startswith("-") else "run"
    if command in ("solve-worker", "handover-worker"):
        from local_optimizer import workflow

        parser = argparse.ArgumentParser(
            description="One local calculation requested by the parent."
        )
        parser.add_argument("workspace", type=Path)
        parser.add_argument(
            "--stage" if command == "solve-worker" else "--directory", required=True
        )
        options = parser.parse_args(args)
        if command == "solve-worker":
            workflow.solve_worker(options.workspace, options.stage)
        else:
            workflow.handover_worker(options.workspace, options.directory)
        return 0
    if command == "plot":
        from local_optimizer import report

        parser = argparse.ArgumentParser(description="Draw the completed workflow's saved values.")
        parser.add_argument("workspace", type=Path)
        parser.add_argument("--output-dir", type=Path)
        options = parser.parse_args(args)
        report.plot(options.workspace, output=options.output_dir)
        return 0
    if command == "plot-saved":
        from local_optimizer import report

        parser = argparse.ArgumentParser(description="Plot a saved schedule; no new flight.")
        parser.add_argument("--result", type=Path, required=True)
        parser.add_argument("--output-dir", type=Path, required=True)
        options = parser.parse_args(args)
        report.draw_saved_result(report.load_saved_result(options.result), options.output_dir)
        return 0

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-name", default="mercury_01")
    parser.add_argument("--data-dir", type=Path)
    options = parser.parse_args(args)
    if command not in ("configure", "check", "run"):
        parser.error("Use configure, check, run, plot or plot-saved")
    from local_optimizer import case, console

    try:
        if command in ("check", "run"):
            print("Checking Mercury inputs and the Python environment...", flush=True)
        prepared = case.configure(WORKFLOW_ROOT.parent, options.run_name, options.data_dir)
        if command == "configure":
            print(
                json.dumps(
                    {"status": "CONFIGURED_NO_SCIENCE", "settings": str(prepared["settings"])}
                )
            )
            return 0
        result = case.dispatch(prepared, command)
        if command == "check":
            print(json.dumps(result, indent=2))
            return 0
        console.print_entry(result)
        return 0 if result["workflow_complete"] else 1
    except (Exception, KeyboardInterrupt) as error:
        print("STOP: " + type(error).__name__ + ": " + str(error), file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
