"""Read one Mercury case and prepare a shared Mac/Windows run."""

from dataclasses import dataclass
import importlib.metadata
import json
import os
from pathlib import Path
import shutil
import struct
import subprocess
import sys
import tempfile
import time
from . import bridge, controls, objective, records, workflow

ROOT = Path(__file__).resolve().parents[2]
SCHEMA = "mercury_shared_case_v1"
RECEIPT = "intrance_native_capabilities_receipt_v1"
require = records.require
CAPABILITIES = {
    "schema": "intrance_unified_v1",
    "modes": ["search", "export", "schedule"],
    "profiles": ["mercury_body_rendezvous_v1", "mercury_body_rendezvous_precision_v3"],
    "observer_schema": 1,
    "native_build": "x86_pe32",
    "pointer_bits": 32,
    "fresh_process_per_stage": True,
}
SETTINGS = {
    "native_sha256",
    "native_executable",
    "data_directory",
    "mission_profile",
    "capabilities_file",
    "seed",
    "plot_python",
    "run_directory",
    "reference_commands",
    "schema",
    "workflow",
    "mode",
    "evaluations",
}


def resolve(value, base, *, exists=True):
    require(
        isinstance(value, str) and bool(value.strip()) and ("\x00" not in value),
        "Expected a file path",
    )
    path = Path(value)
    if not path.is_absolute():
        path = Path(base) / path
    path = path.resolve(strict=exists)
    return path


def runtime():
    versions = {name: importlib.metadata.version(name) for name in ("casadi", "numpy")}
    require(sys.version_info[:2] == (3, 12), "Use the reviewed Python 3.12 environment")
    require(
        versions == {"casadi": "3.7.2", "numpy": "2.5.3"},
        "Install the pinned CasADi/NumPy versions",
    )
    return {
        "python": sys.version.split()[0],
        "pointer_bits": struct.calcsize("P") * 8,
        "packages": versions,
        "check": "package metadata; no NLP or flight run",
    }


def interpreter(value, base):
    path = Path(sys.executable) if value is None else Path(value)
    if not path.is_absolute():
        path = Path(base) / path
    path = Path(os.path.abspath(path))
    resolve(str(path), base)
    return path


def _valid_capabilities(value):
    return type(value) is dict and json.dumps(value, sort_keys=True, allow_nan=False) == json.dumps(
        CAPABILITIES, sort_keys=True, allow_nan=False
    )


@dataclass(frozen=True)
class BuildReceipt:
    executable_sha256: str
    capabilities: dict[str, object]


def _read_build_receipt(path) -> BuildReceipt:
    receipt = records.load(path)
    require(
        type(receipt) is dict
        and set(receipt) == {"schema", "executable_sha256", "capabilities"}
        and (receipt["schema"] == RECEIPT),
        "Wrong build capability receipt",
    )
    digest = receipt["executable_sha256"]
    require(
        isinstance(digest, str) and bridge.HEX.fullmatch(digest), "Wrong build executable SHA-256"
    )
    require(
        _valid_capabilities(receipt["capabilities"]),
        "Native receipt does not support this unified workflow",
    )
    return BuildReceipt(digest, receipt["capabilities"])


def plot_check(interpreter):
    script = (
        "import json,numpy,matplotlib;matplotlib.use('Agg');from "
        "matplotlib.figure import Figure;from "
        "matplotlib.backends.backend_agg import FigureCanvasAgg;f=Fig"
        "ure(figsize=(1,1));c=FigureCanvasAgg(f);f.add_subplot().plot"
        "([0,1],[0,1]);c.draw();print(json.dumps({'numpy':numpy.__ver"
        "sion__,'matplotlib':matplotlib.__version__,'agg':True}))"
    )
    with tempfile.TemporaryDirectory(prefix="intrance_entry_plot_") as cache:
        environment = dict(
            os.environ,
            MPLBACKEND="Agg",
            MPLCONFIGDIR=cache,
            OMP_NUM_THREADS="1",
            OPENBLAS_NUM_THREADS="1",
            MKL_NUM_THREADS="1",
        )
        for key in ("PYTHONHOME", "PYTHONPATH", "PYTHONUSERBASE", "PYTHONNOUSERSITE"):
            environment.pop(key, None)
        done = subprocess.run(
            [str(interpreter), "-B", "-c", script],
            cwd=cache,
            env=environment,
            capture_output=True,
            timeout=30,
            check=False,
        )
    require(
        done.returncode == 0,
        "Plot environment check failed: " + done.stderr.decode(errors="replace")[-1000:],
    )
    result = json.loads(done.stdout)
    require(
        type(result) is dict
        and set(result) == {"numpy", "matplotlib", "agg"}
        and (result["agg"] is True),
        "Incomplete plot check",
    )
    return result


def _probe_capabilities(executable, digest, receipt):
    require(receipt is not None, "Run CHECK_NATIVE.cmd once in Windows for this folder")
    require(
        receipt.executable_sha256 == digest, "Native setup record belongs to another executable"
    )
    return {"schema": RECEIPT, "executable_sha256": digest, "capabilities": receipt.capabilities}


def check(input_file, settings_file):
    settings_file = resolve(str(settings_file), Path.cwd())
    input_file = resolve(str(input_file), Path.cwd())
    settings = records.load(settings_file)
    require(type(settings) is dict, "Local settings must be a JSON object")
    require(
        set(settings) == SETTINGS and settings.get("schema") == SCHEMA,
        "Wrong local-settings fields/schema",
    )
    require(settings["mode"] == "shared", "This workflow uses the shared Mac/Windows route")
    require(settings["workflow"] == "local_refinement", "Use the Mercury local-refinement route")
    require(input_file.name == "mercury.inp", "This entry supports the fixed mercury.inp case")
    for name in bridge.INPUTS:
        source = input_file.parent / name
        require(
            source.is_file() and (not source.is_symlink()),
            "Missing or linked native input: " + name,
        )
    bridge._fixed_case(input_file.parent)
    bridge._native_input(input_file.read_bytes(), settings["seed"], True)
    base = settings_file.parent
    paths = {
        name: resolve(settings[name], base)
        for name in ("native_executable", "data_directory", "mission_profile")
    }
    paths["plot_python"] = interpreter(settings["plot_python"], base)
    paths["run_directory"] = resolve(settings["run_directory"], base, exists=False)
    paths["reference_commands"] = resolve(settings["reference_commands"], base)
    require(
        controls.build_policy(paths["reference_commands"]) == controls.frozen_reference_policy(),
        "Reference commands do not match the command-change limit",
    )
    require(paths["run_directory"].parent.is_dir(), "Create the parent output directory first")
    require(
        not paths["run_directory"].is_relative_to(ROOT),
        "Keep run directories outside the source checkout",
    )
    require(
        paths["native_executable"].is_file()
        and paths["plot_python"].is_file()
        and paths["data_directory"].is_dir(),
        "Wrong executable/data/plot path type",
    )
    receipt_path = (
        None
        if settings["capabilities_file"] is None
        else resolve(settings["capabilities_file"], base)
    )
    build_receipt = None if receipt_path is None else _read_build_receipt(receipt_path)
    digest = settings["native_sha256"]
    if digest is None:
        require(
            build_receipt is not None, "A null native_sha256 requires the build capability receipt"
        )
        digest = build_receipt.executable_sha256
    else:
        require(
            isinstance(digest, str) and bridge.HEX.fullmatch(digest),
            "Unified native executable SHA-256 differs",
        )
    require(
        records.sha(paths["native_executable"]) == digest,
        "Unified native executable SHA-256 differs",
    )
    for name in bridge.DATA_FILES:
        source = paths["data_directory"] / name
        require(
            source.is_file() and (not source.is_symlink()), "Missing or linked native data: " + name
        )
    mission = records.load(paths["mission_profile"])
    require(type(mission) is dict, "Mission profile file must contain a JSON object")
    profile = workflow.profile(settings["seed"], settings["evaluations"], mission)
    receipt = _probe_capabilities(paths["native_executable"], digest, build_receipt)
    result = {
        "schema": SCHEMA,
        "settings": settings,
        "settings_file": str(settings_file),
        "workflow_name": settings["workflow"],
        "settings_sha256": records.sha(settings_file),
        "input_file": str(input_file),
        "input_hashes": {name: records.sha(input_file.parent / name) for name in bridge.INPUTS},
        "paths": {key: str(value) for key, value in paths.items()},
        "native_sha256": digest,
        "native_capabilities": receipt,
        "mission_profile": profile["mission_profile"],
        "runtime": runtime(),
        "plot_check": plot_check(paths["plot_python"]),
        "maximum_solver_calls": profile["maximum_solver_calls"],
        "scientific_calls": 0,
    }
    result["reference_commands_sha256"] = records.sha(paths["reference_commands"])
    return result


def _entry_pins(work, sidecar):
    files = [work / "entry.json", *sidecar.rglob("*")]
    return records.fingerprint([p for p in files if p.is_file()] + [ROOT / "run.py"])


def prepare(checked):
    work = Path(checked["paths"]["run_directory"])
    sidecar = work.with_name(work.name + ".inputs")
    require(
        not work.exists() and (not sidecar.exists()),
        "Choose a new run name; old outputs and .inputs are preserved",
    )
    sidecar.mkdir()
    inputs = sidecar / "native_inputs"
    inputs.mkdir()
    source = Path(checked["input_file"]).parent
    for name in bridge.INPUTS:
        require(
            records.sha(source / name) == checked["input_hashes"][name],
            "Native input changed during checking",
        )
        shutil.copy2(source / name, inputs / name)
    require(
        records.sha(Path(checked["settings_file"])) == checked["settings_sha256"],
        "Local settings changed during checking",
    )
    shutil.copy2(checked["settings_file"], sidecar / "settings.json")
    records.save(sidecar / "capabilities.json", checked["native_capabilities"])
    records.save(sidecar / "mission.json", checked["mission_profile"])
    paths = checked["paths"]
    source_reference = Path(paths["reference_commands"])
    require(
        records.sha(source_reference) == checked["reference_commands_sha256"],
        "Reference commands changed after checking",
    )
    shutil.copy2(source_reference, sidecar / "reference_commands.csv")
    require(
        records.sha(sidecar / "reference_commands.csv") == checked["reference_commands_sha256"],
        "Reference command copy changed during preparation",
    )
    installation = {
        "reference_inputs": str(inputs),
        "data_directory": paths["data_directory"],
        "plot_python": paths["plot_python"],
    }
    for role in ("search", "export", "precision"):
        installation[role + "_executable"] = paths["native_executable"]
        installation[role + "_sha256"] = checked["native_sha256"]
    installation_path = sidecar / "installation.json"
    records.save(installation_path, installation)
    settings = checked["settings"]
    options = {
        "seed": settings["seed"],
        "evaluation_limit": settings["evaluations"],
        "mission_profile": checked["mission_profile"],
    }
    prepared = workflow.prepare(
        installation_path, work, sidecar / "reference_commands.csv", **options
    )
    records.save(work / "entry.json", checked)
    binding = records.load(work / "_audit/inputs.json")
    binding.update(_entry_pins(work, sidecar))
    temporary = work / "_audit/entry_inputs.tmp"
    records.save(temporary, binding)
    temporary.replace(work / "_audit/inputs.json")
    return {
        **prepared,
        "entry_mode": settings["mode"],
        "workflow_name": checked["workflow_name"],
        "entry_inputs": str(sidecar),
    }


def _prepared(checked):
    work = Path(checked["paths"]["run_directory"])
    record = records.load(work / "entry.json")
    require(
        record.get("workflow_name") == "local_refinement",
        "Prepared workflow differs; use a new run folder",
    )
    require(
        record["settings_sha256"] == checked["settings_sha256"]
        and record["input_hashes"] == checked["input_hashes"]
        and (record["paths"] == checked["paths"])
        and (record["native_sha256"] == checked["native_sha256"])
        and (record["native_capabilities"] == checked["native_capabilities"]),
        "Prepared entry settings, inputs or paths changed",
    )
    recorded_profile = objective.validate_profile(record["mission_profile"])
    current_profile = checked["mission_profile"]
    require(
        recorded_profile == current_profile, "Prepared mission profile changed; choose a new run"
    )
    require(
        record["runtime"] == checked["runtime"],
        "Prepared numerical runtime changed; choose a new run",
    )
    require(
        record.get("reference_commands_sha256") == checked["reference_commands_sha256"],
        "Prepared reference commands changed; choose a new run",
    )
    config = workflow.checked(work)
    require(
        config["mission_profile"] == current_profile,
        "Prepared workflow mission profile differs from the requested profile",
    )
    require(
        config["installation"]["search_executable"]
        == config["installation"]["export_executable"]
        == config["installation"]["precision_executable"]
        == checked["paths"]["native_executable"],
        "Prepared entry no longer uses one unified executable",
    )
    require(
        records.sha(work / "_audit/reference_commands.csv") == checked["reference_commands_sha256"],
        "Prepared guard-reference copy differs from the requested file",
    )
    return work


def run(checked):
    began = time.monotonic()
    work = Path(checked["paths"]["run_directory"])
    new = not work.exists()
    if new:
        prepare(checked)
    _prepared(checked)
    require(
        not any(
            (
                (work / name).exists()
                for name in ("attempt.json", "entry_attempt.json", "entry_result.json")
            )
        ),
        "This run was already attempted; choose a new name",
    )
    records.save(
        work / "entry_attempt.json",
        {
            "schema": SCHEMA,
            "mode": "shared",
            "workflow_name": "local_refinement",
            "started_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        },
    )
    outcome, error = (None, None)
    try:
        print("Open START_WINDOWS_WORKER.cmd from this same folder in Windows.", flush=True)
        print(work / "native/RUN_NATIVE_WORKER.cmd", flush=True)
        outcome = workflow.run(work, started=began if new else None)
    except (Exception, KeyboardInterrupt) as exc:
        error = type(exc).__name__ + ": " + str(exc)
    result = {
        "schema": "intrance_local_entry_result_v2",
        "mode": "shared",
        "workflow_name": "local_refinement",
        "workflow_complete": bool(
            outcome and outcome.get("workflow_complete") is True and (error is None)
        ),
        "error": error,
        "worker_cleanup": None,
        "entry_elapsed_seconds": time.monotonic() - began,
        "prepared_in_this_call": new,
        "workflow_result": outcome,
    }
    records.save(work / "entry_result.json", result)
    return result


import re

DATA_FILES = bridge.DATA_FILES
NATIVE_INPUTS = bridge.INPUTS
CASE = Path("case")
RESERVED_NAMES = {
    "CON",
    "PRN",
    "AUX",
    "NUL",
    *(f"COM{i}" for i in range(1, 10)),
    *(f"LPT{i}" for i in range(1, 10)),
}


def run_name(value):
    if not re.fullmatch("[A-Za-z0-9][A-Za-z0-9_-]{0,63}", value) or value.upper() in RESERVED_NAMES:
        raise ValueError(
            "Use a run name with letters, numbers, underscores or hyphens (1..64 characters)"
        )
    return value


def portable_path(path, base):
    try:
        return Path(os.path.relpath(path, base)).as_posix()
    except ValueError:
        return str(path)


def configure(bundle, name, data_directory=None):
    """Create new-run settings only; no scientific calculation."""
    mode = "shared"
    name = run_name(name)
    bundle = Path(bundle).resolve(strict=True)
    workflow_root = bundle / "workflow"
    settings_dir = bundle / "settings"
    output = bundle / "runs" / mode / name
    if not workflow_root.is_dir() or not (workflow_root / "run.py").is_file():
        raise ValueError("Use the extracted preview bundle with its workflow/ directory")
    if output.exists() or output.with_name(output.name + ".inputs").exists():
        raise ValueError("This run directory already exists; choose a new --run-name")
    data = Path(data_directory) if data_directory is not None else bundle / "native_data"
    if not data.is_absolute():
        data = bundle / data
    data = data.resolve(strict=True)
    receipt = bundle / "_internal/runtime/InTrance.capabilities.json"
    if not receipt.is_file():
        raise ValueError(
            "Run CHECK_NATIVE.cmd once in Windows to check this installed EXE before setup"
        )
    required = [
        workflow_root / "runtime/InTrance.exe",
        bundle / "_internal/runtime/InTrance.capabilities.json",
        workflow_root / "case/arrival_limits.json",
        workflow_root / "case/reference_commands.csv",
    ]
    required += [workflow_root / CASE / item for item in NATIVE_INPUTS]
    required += [data / item for item in DATA_FILES]
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise ValueError("Missing release/input files: " + ", ".join(missing))
    value = {
        "schema": SCHEMA,
        "mode": mode,
        "workflow": "local_refinement",
        "native_executable": portable_path(workflow_root / "runtime/InTrance.exe", settings_dir),
        "native_sha256": None,
        "capabilities_file": portable_path(
            bundle / "_internal/runtime/InTrance.capabilities.json", settings_dir
        ),
        "data_directory": portable_path(data, settings_dir),
        "run_directory": portable_path(output, settings_dir),
        "mission_profile": portable_path(workflow_root / "case/arrival_limits.json", settings_dir),
        "reference_commands": portable_path(
            workflow_root / "case/reference_commands.csv", settings_dir
        ),
        "plot_python": None,
        "seed": 1,
        "evaluations": 200,
    }
    destination = settings_dir / f"{mode}-{name}.local.private.json"
    if destination.exists():
        if json.loads(destination.read_text(encoding="utf-8")) != value:
            raise ValueError(
                "These settings already exist with different paths; choose a new run name"
            )
    else:
        settings_dir.mkdir(parents=True, exist_ok=True)
        with destination.open("x", encoding="utf-8") as stream:
            json.dump(value, stream, indent=2, allow_nan=False)
            stream.write("\n")
    output.parent.mkdir(parents=True, exist_ok=True)
    return {
        "app": workflow_root,
        "settings": destination,
        "input": workflow_root / CASE / "mercury.inp",
        "output": output,
        "mode": mode,
        "run_name": name,
    }


def dispatch(prepared, command):
    checked = check(prepared["input"], prepared["settings"])
    if command == "check":
        return {
            "status": "CHECKED_NO_SCIENCE",
            "mode": "shared",
            "run_directory": str(prepared["output"]),
            "runtime": checked["runtime"],
            "plot_check": checked["plot_check"],
        }
    pointer = prepared["settings"].parent / "shared_worker.txt"
    temporary = pointer.with_suffix(".tmp")
    temporary.write_text(prepared["run_name"] + "\n", encoding="ascii")
    os.replace(temporary, pointer)
    return run(checked)
