"""Prepare four native stages and exchange their requests/results with Windows."""

import math
from pathlib import Path
import re
import shutil
import time
import uuid

from . import native_results
from .native_results import BANK_OUTPUTS, search_contract
from .records import require, sha, load, atomic_new, input_file, hash_bytes


INPUTS = ("mercury.inp", "mercury.sim", "mercury.scp", "mercury.ncp", "mercury.eap")


DATA_FILES = ("bin1600_600.405", "ELEMENTS.UNNUM", "ELEMENTS.NUMBR", "ELEMENTS.COMET")


SCHEDULE_FILES = ("bridge_input.csv", "bridge_initial_input.csv", "bridge_schedule_limits.csv")


V1 = "mercury_body_rendezvous_v1"


V3 = "mercury_body_rendezvous_precision_v3"


HEX = re.compile(r"[0-9a-f]{64}\Z")


NONCE = re.compile(r"[0-9a-f]{32}\Z")


BASE_OUTPUTS = (
    "bridge_initial.csv",
    "bridge_states.csv",
    "bridge_commands.csv",
    "bridge_rhs.csv",
    "mercury.csv",
    "mercury.ctr",
    "mercury.dat",
    "mercury.gesop.txt",
    "mercury.wrl",
    "mercury.rep",
    "stdout.txt",
    "stderr.txt",
    "native_exit_code.txt",
    "guest_identity_verified.txt",
)


PRECISE_OUTPUTS = ("bridge_endpoint.csv", "bridge_target_endpoint.csv", "bridge_intervals.csv")


STAGES = ("cold_search", "export", "precise_baseline", "precise_candidate")


SESSION_SECONDS, CHILD_SECONDS = 7200, 90


SCHEMA = 5


KIND = "fresh_enc_workflow_transport_v3"


REFERENCE_INPUTS = Path(__file__).resolve().parents[2] / "case"


SEARCH_OUTPUTS = (
    "mercury.eac",
    "mercury.csv",
    "mercury.ctr",
    "mercury.dat",
    "mercury.gesop.txt",
    "mercury.wrl",
    "mercury.rep",
    "stdout.txt",
    "stderr.txt",
    "native_exit_code.txt",
    "guest_identity_verified.txt",
)


WORKER_FILE = Path(__file__).with_name("native_worker.ps1")


def _checked_native_input(data):
    """Only the six supplied mission filenames and four sibling data paths."""
    lines = data.decode("utf-8").splitlines()
    parsed = {}
    for line in lines:
        clean = line.split("//", 1)[0].strip()
        if not clean:
            continue
        require("=" in clean, "Malformed INP line")
        key, value = [p.strip() for p in clean.split("=", 1)]
        require(key not in parsed, "Duplicate INP field")
        parsed[key] = value
    paths = {
        "DE405_DATA_FILE": "../data/bin1600_600.405",
        "ASTEROID_DATA_FILE1": "../data/ELEMENTS.UNNUM",
        "ASTEROID_DATA_FILE2": "../data/ELEMENTS.NUMBR",
        "COMET_DATA_FILE": "../data/ELEMENTS.COMET",
        "SIM_PARAM_FILE": "mercury.sim",
        "SC_PARAM_FILE": "mercury.scp",
        "NC_PARAM_FILE": "mercury.ncp",
        "EA_PARAM_FILE": "mercury.eap",
        "SIM_DATA_FILE": "mercury.csv",
        "TRAJ_DATA_FILE": "mercury.dat",
        "GESOP_FILE": "mercury.gesop.txt",
        "VRML_FILE": "mercury.wrl",
        "CTRL_FILE": "mercury.ctr",
        "BEST_CHROM_FILE": "mercury.eac",
        "REPORT_FILE": "mercury.rep",
    }
    require(
        set(parsed) == set(paths) | {"FILE_TYPE", "RANDOM_SEED", "COMMAND", "COLDSTART"}
        and all(parsed[k] == v for k, v in paths.items()),
        "Unexpected native input fields/paths",
    )
    require(
        parsed["COMMAND"] in ("optimize NC", "evaluate NC")
        and parsed["RANDOM_SEED"] == "1"
        and parsed["COLDSTART"] == "yes",
        "Only this explicit saved-controller evaluation is permitted",
    )
    lines = [
        "COMMAND=evaluate NC" if line.split("=", 1)[0].strip() == "COMMAND" else line
        for line in lines
    ]
    return ("\r\n".join(lines) + "\r\n").encode("utf-8")


def _simulation_settings(data, precise):
    for key in (b"MAX_RELATIVE_ERROR", b"MAX_ABSOLUTE_ERROR"):
        pattern = rb"(?m)^(" + key + rb"=)1e-6(\r?)$"
        require(len(re.findall(pattern, data)) == 1, "Expected original matched 1e-6 tolerances")
        if precise:
            data = re.sub(pattern, lambda m: m[1] + b"1e-14" + m[2], data)
    return data


def _fixed_case(reference):
    """This release varies only seed/cutoff, not the physical/scoring contract."""
    reference = Path(reference)
    for name in INPUTS:
        require(
            (reference / name).read_bytes() == (REFERENCE_INPUTS / name).read_bytes(),
            "Unsupported mission/ENC input change: " + name,
        )

    def fields(path):
        items = {}
        for raw in path.read_text().splitlines():
            line = raw.split("//", 1)[0].strip()
            if line:
                key, value = (part.strip() for part in line.split("=", 1))
                require(key not in items, "Repeated native setting")
                items[key] = value
        return items

    simulation, evolution = fields(reference / "mercury.sim"), fields(reference / "mercury.eap")
    exact = {
        "TARGET_BODY_NAME": "mercury",
        "TARGET_STATE": "body rendezvous",
        "INITIAL_BODY_NAME": "earth",
        "INITIAL_STATE": "body",
        "OPTIMIZATION_GOAL": "minimum transfer time",
        "INTEGRATOR": "RKF54",
    }
    numbers = {
        "SIM_START_TIME_MIN": 52654.5,
        "SIM_START_TIME_MAX": 52654.5,
        "INTEGRATION_INTERVAL": 600.0,
        "INTEGRATION_STEPS": 600.0,
        "TARGET_DIST_MAX": 100000.0,
        "TARGET_RELVEL_MAX": 0.1,
        "INITIAL_VINF_MIN": 0.0,
        "INITIAL_VINF_MAX": 0.0,
        "MIN_SOLAR_DISTANCE": 0.2,
    }
    require(
        all(simulation.get(key) == value for key, value in exact.items())
        and all(float(simulation.get(key, "nan")) == value for key, value in numbers.items())
        and evolution.get("FITNESS_FUNCTION_TYPE") == "J_and",
        "Fixed Mercury model/native fitness changed",
    )


def _native_input(data, seed, search):
    # Validate every path/field against the fixed case, then alter only command/seed.
    evaluated = _checked_native_input(data)
    command = b"optimize NC" if search else b"evaluate NC"
    return evaluated.replace(b"COMMAND=evaluate NC", b"COMMAND=" + command).replace(
        b"RANDOM_SEED=1\r\n", (f"RANDOM_SEED={seed}\r\n").encode("ascii")
    )


def _observer(config):
    # The final immutable snapshot is always written; no redundant intermediate bank.
    return (
        (
            "schema=1\nseed={seed}\neval_limit={evaluation_limit}\ncapture=false\n"
            "checkpoints={evaluation_limit}\noutput_directory=enc_observer\n"
        )
        .format(**config["search"])
        .encode("ascii")
    )


def _fixture_bytes(root, config, stage_id):
    precise = stage_id.startswith("precise_")
    files = {n: input_file(root, "assets/reference/" + n).read_bytes() for n in INPUTS}
    files["mercury.inp"] = _native_input(
        files["mercury.inp"], config["search"]["seed"], stage_id == "cold_search"
    )
    files["mercury.sim"] = _simulation_settings(files["mercury.sim"], precise)
    if stage_id == "cold_search":
        files["enc_observer.cfg"] = _observer(config)
    return files


def _worker(config_sha):
    return (
        WORKER_FILE.read_text(encoding="ascii")
        .replace("@CONFIG_SHA@", config_sha)
        .replace("\n", "\r\n")
        .encode("ascii")
    )


def configure(
    native_root,
    reference,
    search_executable,
    export_executable,
    precision_executable,
    data_dir,
    search_sha,
    export_sha,
    precision_sha,
    *,
    seed,
    evaluation_limit,
    search_seconds=720,
    terminal_policy="exact_completed_count",
):
    """Prepare one repeatable attempt, without launching and without saved controllers."""
    root, reference, data_dir = map(Path, (native_root, reference, data_dir))
    require(not root.exists(), "Preserve existing native workspace")
    search = search_contract(seed, evaluation_limit, search_seconds, terminal_policy)
    require(
        not (reference / "mercury.eac").exists(),
        "Fresh inputs must not contain a historical controller",
    )
    # The case folder also holds local settings; only INPUTS below enter the native run.
    _fixed_case(reference)
    _native_input((reference / "mercury.inp").read_bytes(), seed, True)
    _simulation_settings((reference / "mercury.sim").read_bytes(), True)
    eap = (reference / "mercury.eap").read_text()
    for key, expected in (("POPULATION_SIZE", "40"), ("SEARCH_SCAN_EPOCHS", "5")):
        found = re.findall(r"(?m)^" + key + r"=([^\r\n]+)", eap)
        require(
            [v.split("//", 1)[0].strip() for v in found] == [expected],
            "Observer requires original population/scan settings",
        )
    source = {
        **{"assets/reference/" + n: reference / n for n in INPUTS},
        **{"assets/data/" + n: data_dir / n for n in DATA_FILES},
    }
    for name, path, expected in (
        ("search", search_executable, search_sha),
        ("export", export_executable, export_sha),
        ("precision", precision_executable, precision_sha),
    ):
        path = Path(path)
        require(
            isinstance(expected, str)
            and HEX.fullmatch(expected)
            and path.is_file()
            and sha(path) == expected,
            "Executable identity mismatch: " + name,
        )
        source["assets/" + name + ".exe"] = path
    require(
        all(p.is_file() and not p.is_symlink() for p in source.values()),
        "Missing/linked source asset",
    )
    before = {str(p.resolve()): sha(p) for p in source.values()}
    root.mkdir(parents=True)
    for relative, path in source.items():
        to = root / relative
        to.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(path, to)
    for name in (
        "stages",
        "_transport/requests",
        "_transport/accepted",
        "_transport/client_claims",
    ):
        (root / name).mkdir(parents=True)
    atomic_new(
        root / "source.json",
        {
            "schema": 1,
            "original_assets": before,
            "controller_supplied": False,
            "build_performed": False,
        },
    )
    config = {
        "schema": SCHEMA,
        "kind": KIND,
        "run_id": uuid.uuid4().hex,
        "stage_ids": list(STAGES),
        "session_seconds": SESSION_SECONDS,
        "child_seconds": CHILD_SECONDS,
        "max_jobs": 4,
        "search": search,
        "source_sha256": sha(root / "source.json"),
        "executables": {"search": search_sha, "export": export_sha, "precision": precision_sha},
        "assets": {n: sha(root / n) for n in sorted(source)},
    }
    config["fixtures"] = {
        stage: {n: hash_bytes(data) for n, data in _fixture_bytes(root, config, stage).items()}
        for stage in STAGES
    }
    atomic_new(root / "transport.json", config)
    (root / "RUN_NATIVE_WORKER.ps1").write_bytes(_worker(sha(root / "transport.json")))
    (root / "RUN_NATIVE_WORKER.cmd").write_bytes(
        LAUNCHER.replace("@WORKER_SHA@", sha(root / "RUN_NATIVE_WORKER.ps1"))
        .replace("\n", "\r\n")
        .encode("ascii")
    )
    require(
        all(sha(Path(p)) == digest for p, digest in before.items()),
        "Original changed during preparation",
    )
    return verify(root)


def _metadata(native_root):
    """Validate the small transport contract, not native payload contents.

    Polling is observational. Full payload hashes remain mandatory at dispatch,
    receipt inspection and finalization; a progress status is not acceptance.
    """
    root = Path(native_root)
    c = load(input_file(root, "transport.json"))
    require(
        set(c)
        == {
            "schema",
            "kind",
            "run_id",
            "stage_ids",
            "session_seconds",
            "child_seconds",
            "max_jobs",
            "search",
            "source_sha256",
            "executables",
            "assets",
            "fixtures",
        }
        and type(c["schema"]) is int
        and c["schema"] == SCHEMA
        and c["kind"] == KIND
        and isinstance(c["run_id"], str)
        and NONCE.fullmatch(c["run_id"])
        and c["stage_ids"] == list(STAGES)
        and type(c["session_seconds"]) is int
        and c["session_seconds"] == SESSION_SECONDS
        and type(c["child_seconds"]) is int
        and c["child_seconds"] == CHILD_SECONDS
        and type(c["max_jobs"]) is int
        and c["max_jobs"] == 4,
        "Wrong fresh transport contract",
    )
    s = c["search"]
    require(
        s
        == search_contract(s["seed"], s["evaluation_limit"], s["child_seconds"], s["stop_policy"]),
        "Search contract changed",
    )
    require(set(c["executables"]) == {"search", "export", "precision"}, "Wrong executable roles")
    names = (
        {"assets/reference/" + n for n in INPUTS}
        | {"assets/data/" + n for n in DATA_FILES}
        | {"assets/" + role + ".exe" for role in c["executables"]}
    )
    require(
        set(c["assets"]) == names
        and all(isinstance(d, str) and HEX.fullmatch(d) for d in c["assets"].values()),
        "Wrong frozen asset bindings",
    )
    require(
        all(
            c["assets"]["assets/" + role + ".exe"] == digest
            for role, digest in c["executables"].items()
        ),
        "Executable binding differs",
    )
    require(
        sha(input_file(root, "source.json")) == c["source_sha256"],
        "Original asset provenance changed",
    )
    require(set(c["fixtures"]) == set(STAGES), "Wrong fixture stage set")
    for stage, fixture in c["fixtures"].items():
        expected = set(INPUTS) | ({"enc_observer.cfg"} if stage == "cold_search" else set())
        require(
            isinstance(fixture, dict)
            and set(fixture) == expected
            and all(isinstance(d, str) and HEX.fullmatch(d) for d in fixture.values()),
            "Wrong fixture bindings",
        )
    worker = _worker(sha(root / "transport.json"))
    require(input_file(root, "RUN_NATIVE_WORKER.ps1").read_bytes() == worker, "Worker changed")
    launcher = (
        LAUNCHER.replace("@WORKER_SHA@", hash_bytes(worker)).replace("\n", "\r\n").encode("ascii")
    )
    require(input_file(root, "RUN_NATIVE_WORKER.cmd").read_bytes() == launcher, "Launcher changed")
    return c


def verify(native_root):
    """Full immutable payload verification; never replaced by status polling."""
    root = Path(native_root)
    c = _metadata(root)
    _fixed_case(root / "assets/reference")
    require(
        all(sha(input_file(root, n)) == digest for n, digest in c["assets"].items()),
        "Frozen asset changed",
    )
    require(
        not (root / "assets/reference/mercury.eac").exists(),
        "Historical controller added to fresh assets",
    )
    require(
        c["fixtures"]
        == {
            stage: {n: hash_bytes(d) for n, d in _fixture_bytes(root, c, stage).items()}
            for stage in STAGES
        },
        "Fixture contract changed",
    )
    return c


def _source(root, config, stage_id):
    # Export must use this run's completed controller, not an old successful trajectory.
    if stage_id == "cold_search":
        return {
            "schema": 1,
            "run_id": config["run_id"],
            "kind": "fresh_cold_search",
            "controller_input": None,
            "seed": config["search"]["seed"],
            "evaluation_limit": config["search"]["evaluation_limit"],
        }
    search = root / "stages/cold_search"
    accepted = inspect(search)
    return {
        "schema": 1,
        "run_id": config["run_id"],
        "kind": "same_run_observer_final",
        "controller_source": "stages/cold_search/result/enc_observer/final.eac",
        "controller_sha256": sha(search / "result/enc_observer/final.eac"),
        "search_receipt_sha256": sha(search / "_audit/receipt.json"),
        "observer_summary_sha256": sha(search / "result/enc_observer/summary.json"),
        "observer_metadata_sha256": sha(search / "result/enc_observer/final.json"),
        "session_id": accepted["session_id"],
        "seed": config["search"]["seed"],
        "completed_evaluations": accepted["observer"]["completed_evaluations"],
    }


def prepare_stage(native_root, stage_id, schedule=None):
    root = Path(native_root)
    config = verify(root)
    require(stage_id in STAGES, "Unknown fresh stage")
    precise = stage_id.startswith("precise_")
    require((schedule is not None) == precise, "Only precise stages require a schedule")
    stage = root / "stages" / stage_id
    require(not stage.exists(), "Preserve existing stage")
    source = _source(root, config, stage_id)
    fixed = _fixture_bytes(root, config, stage_id)
    if stage_id != "cold_search":
        fixed["mercury.eac"] = input_file(root, source["controller_source"]).read_bytes()
    if precise:
        fixed.update({n: input_file(Path(schedule), n).read_bytes() for n in SCHEDULE_FILES})
    for name in ("fixture", "data", "_audit"):
        (stage / name).mkdir(parents=True, exist_ok=False)
    role = "search" if stage_id == "cold_search" else "precision" if precise else "export"
    shutil.copy2(root / "assets" / (role + ".exe"), stage / "InTrance.exe")
    for n, data in fixed.items():
        (stage / "fixture" / n).write_bytes(data)
    for n in DATA_FILES:
        shutil.copy2(root / "assets/data" / n, stage / "data" / n)
    atomic_new(stage / "_audit/source.json", source)
    contract = _stage_contract(config, stage_id)
    manifest = {
        "schema": SCHEMA,
        "stage_id": stage_id,
        **contract,
        "run_id": config["run_id"],
        "config_sha256": sha(root / "transport.json"),
        "source_sha256": sha(stage / "_audit/source.json"),
        "hashes": {
            p.relative_to(stage).as_posix(): sha(p)
            for folder in ("fixture", "data")
            for p in sorted((stage / folder).iterdir())
        },
    }
    manifest["hashes"]["InTrance.exe"] = sha(stage / "InTrance.exe")
    atomic_new(stage / "_audit/manifest.json", manifest)
    verify_stage(stage)
    return stage


def _stage_contract(config, stage_id):
    search, precise = stage_id == "cold_search", stage_id.startswith("precise_")
    role = "search" if search else "precision" if precise else "export"
    return {
        "mode": "cold_search" if search else "schedule" if precise else "export",
        "profile": "enc_count_observer_v1" if search else V3 if precise else V1,
        "tolerance": "1e-14" if precise else "1e-6",
        "command": "optimize NC" if search else "evaluate NC",
        "executable_sha256": config["executables"][role],
        "child_seconds": config["search"]["child_seconds"] if search else CHILD_SECONDS,
    }


def verify_stage(stage):
    stage = Path(stage)
    root = stage.parent.parent
    config = verify(root)
    require(
        stage.parent.name == "stages" and stage.name in STAGES and not stage.is_symlink(),
        "Wrong stage path",
    )
    m = load(input_file(stage, "_audit/manifest.json"))
    contract = _stage_contract(config, stage.name)
    require(
        set(m)
        == {"schema", "stage_id", "run_id", "config_sha256", "source_sha256", "hashes"}
        | set(contract)
        and type(m["schema"]) is int
        and m["schema"] == SCHEMA
        and m["stage_id"] == stage.name
        and m["run_id"] == config["run_id"]
        and m["config_sha256"] == sha(root / "transport.json")
        and all(type(m[k]) is type(v) and m[k] == v for k, v in contract.items()),
        "Wrong stage contract",
    )
    require(
        sha(input_file(stage, "_audit/source.json")) == m["source_sha256"]
        and load(stage / "_audit/source.json") == _source(root, config, stage.name),
        "Same-run controller lineage differs",
    )
    files = set(config["fixtures"][stage.name])
    if stage.name != "cold_search":
        files.add("mercury.eac")
    if stage.name.startswith("precise_"):
        files.update(SCHEDULE_FILES)
    names = {"InTrance.exe"} | {"fixture/" + n for n in files} | {"data/" + n for n in DATA_FILES}
    require(
        set(m["hashes"]) == names
        and all(sha(input_file(stage, n)) == h for n, h in m["hashes"].items()),
        "Stage payload differs",
    )
    require(
        {p.name for p in (stage / "fixture").iterdir()} == files, "Unexpected native input file"
    )
    require(m["hashes"]["InTrance.exe"] == m["executable_sha256"], "Wrong stage executable")
    require(
        all(m["hashes"]["fixture/" + n] == h for n, h in config["fixtures"][stage.name].items()),
        "Fixed mission input changed",
    )
    require(
        all(m["hashes"]["data/" + n] == config["assets"]["assets/data/" + n] for n in DATA_FILES),
        "Planetary data changed",
    )
    if stage.name != "cold_search":
        require(
            m["hashes"]["fixture/mercury.eac"]
            == sha(root / "stages/cold_search/result/enc_observer/final.eac"),
            "Historical controller substitution",
        )
    return m


def status(native_root):
    """Read progress with metadata checks only; use inspect() for acceptance."""
    root = Path(native_root)
    config = _metadata(root)
    if not (root / "_transport/session.json").exists():
        return {"state": "NOT_STARTED", "next_stage": STAGES[0]}
    session = load(input_file(root, "_transport/session.json"))
    require(
        set(session)
        == {"schema", "session_id", "run_id", "config_sha256", "worker_pid", "session_seconds"}
        and type(session["schema"]) is int
        and session["schema"] == SCHEMA
        and isinstance(session["session_id"], str)
        and NONCE.fullmatch(session["session_id"])
        and session["run_id"] == config["run_id"]
        and session["config_sha256"] == sha(root / "transport.json")
        and type(session["worker_pid"]) is int
        and session["worker_pid"] > 0
        and type(session["session_seconds"]) is int
        and session["session_seconds"] == SESSION_SECONDS,
        "Invalid fresh worker session",
    )
    files = sorted((root / "_transport").glob("status_*.json"))
    if not files:
        return {**session, "state": "STARTING"}
    require(
        all(re.fullmatch(r"status_[0-9]{6}\.json", p.name) for p in files),
        "Invalid status filename",
    )
    record = load(input_file(root, "_transport/" + files[-1].name))
    require(
        set(record)
        == {
            "schema",
            "session_id",
            "state",
            "sequence",
            "challenge",
            "elapsed_seconds",
            "failure",
            "termination",
            "owned_child_cleanup_verified",
        }
        and type(record["schema"]) is int
        and record["schema"] == SCHEMA
        and record["session_id"] == session["session_id"]
        and record["state"] in ("WAITING", "RUNNING", "COMPLETED", "FAILED", "CANCELLED")
        and type(record["sequence"]) is int
        and 0 <= record["sequence"] <= 4
        and isinstance(record["challenge"], str)
        and NONCE.fullmatch(record["challenge"])
        and type(record["elapsed_seconds"]) in (int, float)
        and math.isfinite(record["elapsed_seconds"])
        and record["elapsed_seconds"] >= 0
        and type(record["owned_child_cleanup_verified"]) is bool,
        "Invalid fresh worker status",
    )
    _check_termination(
        record["state"],
        record["termination"],
        record["failure"],
        record["owned_child_cleanup_verified"],
    )
    require(
        record["state"] != "COMPLETED" or record["sequence"] == len(STAGES),
        "Worker completed before all ordered stages",
    )
    require(
        record["state"] not in ("WAITING", "RUNNING") or record["sequence"] < len(STAGES),
        "Worker waits/runs beyond final stage",
    )
    if record["state"] in ("WAITING", "RUNNING", "COMPLETED") or record["sequence"] > 0:
        _check_runtime_preflight(root, session, config)
    return {
        **session,
        **record,
        "next_stage": STAGES[record["sequence"]] if record["sequence"] < 4 else None,
    }


def runtime_preflight(native_root, *, session=None):
    """Inspect the real Windows preflight receipt; no executable is launched here."""
    root = Path(native_root)
    config = verify(root)
    if session is None:
        session = load(input_file(root, "_transport/session.json"))
    return _check_runtime_preflight(root, session, config)


def _check_runtime_preflight(root, session, config):
    # File sizes are inexpensive metadata; no payload hashing occurs in polls.
    record = load(input_file(root, "_transport/preflight.json"))
    fields = {
        "schema",
        "run_id",
        "session_id",
        "powershell_version",
        "os_version",
        "shared_directory",
        "local_directory",
        "asset_bytes",
        "required_local_bytes",
        "available_local_bytes",
        "native_images",
        "syntax_checked",
        "shared_session_published",
        "ready",
    }
    require(
        set(record) == fields
        and type(record["schema"]) is int
        and record["schema"] == SCHEMA
        and record["run_id"] == config["run_id"]
        and record["session_id"] == session["session_id"]
        and record["native_images"]
        == {role: "x86_pe32" for role in ("search", "export", "precision")}
        and all(record[k] is True for k in ("syntax_checked", "shared_session_published", "ready")),
        "Native runtime preflight failed or belongs to another run",
    )
    require(
        all(
            isinstance(record[k], str) and bool(record[k])
            for k in ("powershell_version", "os_version", "shared_directory", "local_directory")
        )
        and re.fullmatch(r"5\.1\.\d+(?:\.\d+)?", record["powershell_version"]),
        "Unexpected Windows PowerShell runtime",
    )
    require(
        all(
            type(record[k]) is int and record[k] > 0
            for k in ("asset_bytes", "required_local_bytes", "available_local_bytes")
        )
        and record["required_local_bytes"] == 4 * record["asset_bytes"] + 64 * 1024 * 1024
        and record["available_local_bytes"] >= record["required_local_bytes"],
        "Invalid native staging-space preflight",
    )
    require(
        record["asset_bytes"] == sum(input_file(root, n).stat().st_size for n in config["assets"]),
        "Native asset-size preflight differs",
    )
    return record


def _check_request(value, session, stage_id, manifest_sha):
    require(
        set(value)
        == {"schema", "session_id", "sequence", "challenge", "stage_id", "manifest_sha256"}
        and type(value["schema"]) is int
        and value["schema"] == SCHEMA
        and value["session_id"] == session
        and type(value["sequence"]) is int
        and value["sequence"] == STAGES.index(stage_id)
        and value["stage_id"] == stage_id
        and value["manifest_sha256"] == manifest_sha
        and isinstance(value["challenge"], str)
        and NONCE.fullmatch(value["challenge"]),
        "Request identity mismatch",
    )


def request(native_root, stage_id):
    root = Path(native_root)
    require(stage_id in STAGES, "Unknown fresh stage")
    stage = root / "stages" / stage_id
    verify_stage(stage)
    current = status(root)
    sequence = STAGES.index(stage_id)
    published = [
        root / f"_transport/{folder}/{sequence}.json" for folder in ("requests", "accepted")
    ]
    require(sum(p.exists() for p in published) <= 1, "Duplicate request state")
    for path in published:
        if path.exists():
            previous = load(input_file(root, path.relative_to(root).as_posix()))
            _check_request(
                previous, current.get("session_id"), stage_id, sha(stage / "_audit/manifest.json")
            )
            return previous
    require(
        current["state"] == "WAITING"
        and current["next_stage"] == stage_id
        and current["elapsed_seconds"] < SESSION_SECONDS,
        "Worker not ready for ordered fresh stage",
    )
    for prior in STAGES[:sequence]:
        inspect(root / "stages" / prior)
    require(
        not (stage / "attempt.claim").exists() and not (stage / "result").exists(),
        "Stage already attempted",
    )
    require(not list((root / "_transport/requests").glob("*.json")), "Pending native request")
    (root / "_transport/client_claims" / str(sequence)).mkdir()
    value = {
        "schema": SCHEMA,
        "session_id": current["session_id"],
        "sequence": sequence,
        "challenge": current["challenge"],
        "stage_id": stage_id,
        "manifest_sha256": sha(stage / "_audit/manifest.json"),
    }
    atomic_new(root / f"_transport/requests/{sequence}.json", value)
    return value


def inspect(stage):
    stage = Path(stage)
    m = verify_stage(stage)
    root = stage.parent.parent
    session = status(root)
    accepted = load(input_file(root, f"_transport/accepted/{STAGES.index(stage.name)}.json"))
    receipt = load(input_file(stage, "_audit/receipt.json"))
    _check_request(
        accepted, session.get("session_id"), stage.name, sha(stage / "_audit/manifest.json")
    )
    fields = {
        "schema",
        "session_id",
        "sequence",
        "challenge",
        "stage_id",
        "manifest_sha256",
        "success",
        "failure",
        "exit_code",
        "timed_out",
        "process_wall_seconds",
        "cpu_seconds",
        "pid",
        "inputs_before",
        "inputs_after",
        "output_hashes",
        "guest_directory",
        "owned_child_cleanup_verified",
        "timing_source",
        "termination",
    }
    require(
        set(receipt) == fields
        and type(receipt["schema"]) is int
        and receipt["schema"] == SCHEMA
        and all(receipt[k] == accepted[k] for k in accepted)
        and receipt["success"] is True
        and receipt["failure"] is None
        and receipt["timed_out"] is False
        and receipt["owned_child_cleanup_verified"] is True
        and type(receipt["exit_code"]) is int
        and receipt["exit_code"] == 0
        and type(receipt["pid"]) is int
        and receipt["pid"] > 0
        and receipt["inputs_before"] == receipt["inputs_after"] == m["hashes"]
        and receipt["timing_source"] == "Stopwatch"
        and receipt["termination"] is None,
        "Native receipt failed or changed",
    )
    require(
        all(
            type(receipt[k]) in (int, float) and math.isfinite(receipt[k]) and receipt[k] >= 0
            for k in ("process_wall_seconds", "cpu_seconds")
        )
        and 0 < receipt["process_wall_seconds"] <= m["child_seconds"],
        "Invalid native timing/deadline",
    )
    result = stage / "result"
    names = set(INPUTS) | set(SEARCH_OUTPUTS)
    if stage.name == "cold_search":
        names |= {"enc_observer.cfg"} | {"enc_observer/" + n for n in BANK_OUTPUTS}
    else:
        names |= set(BASE_OUTPUTS)
        if stage.name.startswith("precise_"):
            names |= set(SCHEDULE_FILES) | set(PRECISE_OUTPUTS)
    actual = {p.relative_to(result).as_posix() for p in result.rglob("*") if p.is_file()}
    require(set(receipt["output_hashes"]) == actual == names, "Wrong native result file set")
    require(
        all(sha(input_file(result, n)) == h for n, h in receipt["output_hashes"].items()),
        "Native result changed",
    )
    for n in set(p.name for p in (stage / "fixture").iterdir()):
        require(
            input_file(result, n).read_bytes() == (stage / "fixture" / n).read_bytes(),
            "Returned input changed",
        )
    require(
        (result / "native_exit_code.txt").read_text().strip() == "0"
        and (result / "guest_identity_verified.txt").stat().st_size > 0,
        "Native identity markers absent",
    )
    stdout, stderr, report = [
        (result / n).read_text() for n in ("stdout.txt", "stderr.txt", "mercury.rep")
    ]
    require(
        stdout.count("InTrance run finished!") == report.count("InTrance run finished!") == 1,
        "Native completion absent/repeated",
    )
    profiles = re.findall(r"NATIVE_BRIDGE_PROFILE_COMPLETE profile=([^\s]+)", stdout)
    require(
        profiles == ([] if stage.name == "cold_search" else [m["profile"]]),
        "Native profile differs",
    )
    require(
        stdout.count("NATIVE_SCHEDULE_REACHED_END") == int(stage.name.startswith("precise_")),
        "Schedule completion differs",
    )
    patterns = (
        r"EXIT_TRACE registry: (?:inserted|remove) key=\d+ pointer=[0-9A-Fa-f]+",
        r"EXIT_TRACE main: (?:before delete intrance|after delete intrance|returning 0)",
        r"EXIT_TRACE CInTrance: before solar system delete [0-9A-Fa-f]+",
        r"EXIT_TRACE CInTrance: solar system deleted, member teardown next",
    )
    require(
        all(
            any(re.fullmatch(pattern, line) for pattern in patterns)
            for line in stderr.splitlines()
            if line.strip()
        )
        and stderr.count("EXIT_TRACE main: returning 0") == 1,
        "Unexpected native stderr",
    )
    observer = (
        native_results.inspect_observer(result, verify(root)["search"])
        if stage.name == "cold_search"
        else None
    )
    summary = None
    if observer is not None:
        summary = {
            "observer": observer,
            "selected_endpoint": native_results.selected_endpoint(result),
            "raw_final_metadata": load(result / "enc_observer/final.json"),
            "native_fitness": observer["best_fitness"],
            "fitness_basis": "observer best controller at final completed boundary",
            "raw_selected_csv": str(result / "mercury.csv"),
            "raw_native_report": str(result / "mercury.rep"),
        }
    return {
        **receipt,
        "native_process_completed": True,
        "mode": m["mode"],
        "profile": m["profile"],
        "observer": observer,
        "search_summary": summary,
        "numerical_acceptance": "pending independent checks",
    }


def wait_ready(native_root, stage_id, timeout=SESSION_SECONDS, poll_seconds=0.5):
    require(
        stage_id in STAGES
        and type(timeout) in (int, float)
        and math.isfinite(timeout)
        and 0 < timeout <= SESSION_SECONDS
        and type(poll_seconds) in (int, float)
        and math.isfinite(poll_seconds)
        and 0 < poll_seconds <= 10,
        "Invalid ready wait",
    )
    started = time.monotonic()
    while time.monotonic() - started < timeout:
        current = status(native_root)
        if current["state"] == "WAITING" and current["next_stage"] == stage_id:
            require(current["elapsed_seconds"] < SESSION_SECONDS, "Worker session deadline reached")
            return current
        require(
            current["state"] not in ("FAILED", "CANCELLED", "COMPLETED")
            and current.get("sequence", 0) <= STAGES.index(stage_id),
            "Worker ended/advanced before stage",
        )
        time.sleep(poll_seconds)
    raise TimeoutError("Fresh worker not ready; no native execution retried")


def wait_result(native_root, stage_id, timeout=None, poll_seconds=0.5, on_progress=None):
    root = Path(native_root)
    require(stage_id in STAGES, "Unknown fresh stage")
    m = verify_stage(root / "stages" / stage_id)
    timeout = m["child_seconds"] + 60 if timeout is None else timeout
    require(
        type(timeout) in (int, float)
        and math.isfinite(timeout)
        and 0 < timeout <= SESSION_SECONDS
        and type(poll_seconds) in (int, float)
        and math.isfinite(poll_seconds)
        and 0 < poll_seconds <= 10,
        "Invalid result wait",
    )
    started = time.monotonic()
    while time.monotonic() - started < timeout:
        if on_progress is not None:
            on_progress()
        stage = root / "stages" / stage_id
        if (stage / "_audit/receipt.json").exists():
            return inspect(stage)
        require(
            status(root)["state"] not in ("FAILED", "CANCELLED", "COMPLETED"),
            "Worker ended without requested result",
        )
        time.sleep(poll_seconds)
    raise TimeoutError("Fresh native result absent; preserve attempt and stop session")


def _stop_reason(reason):
    require(
        isinstance(reason, str) and re.fullmatch(r"[a-z][a-z0-9_]{0,63}", reason),
        "Invalid controlled host stop reason",
    )
    return reason


def _check_termination(state, termination, failure, cleanup_verified):
    if state in ("WAITING", "RUNNING"):
        require(termination is None and failure is None, "Running worker has a termination record")
        return
    require(
        isinstance(termination, dict) and set(termination) == {"kind", "code"},
        "Terminal classification absent",
    )
    _stop_reason(termination["code"])
    expected = {
        "COMPLETED": "completed",
        "CANCELLED": "controlled_host_stop",
        "FAILED": "native_or_transport_error",
    }
    require(
        termination["kind"] == expected[state],
        "Terminal classification conflicts with worker state",
    )
    if state in ("COMPLETED", "CANCELLED"):
        require(
            failure is None and cleanup_verified is True,
            "Controlled terminal state has an unresolved failure",
        )
    else:
        require(isinstance(failure, str) and bool(failure), "Failed worker has no diagnostic")


def finish(native_root, *, reason="host_requested_stop"):
    """Ask the owned worker to stop; a bounded host decision is not a native fault."""
    reason = _stop_reason(reason)
    root = Path(native_root)
    current = status(root)
    if current["state"] in ("NOT_STARTED", "COMPLETED", "FAILED", "CANCELLED"):
        return current
    value = {"schema": SCHEMA, "session_id": current["session_id"], "reason": reason}
    path = root / "_transport/stop.json"
    if path.exists():
        require(load(path) == value, "Stop request belongs to another session")
    else:
        atomic_new(path, value)
    return {"state": "STOP_REQUESTED", "session_id": current["session_id"], "reason": reason}


def stop_and_wait(
    native_root, *, request_stop=True, timeout=10.0, poll_seconds=0.2, reason="host_requested_stop"
):
    reason = _stop_reason(reason)
    require(
        type(request_stop) is bool
        and type(timeout) in (int, float)
        and math.isfinite(timeout)
        and 0 < timeout <= 30
        and type(poll_seconds) in (int, float)
        and math.isfinite(poll_seconds)
        and 0 < poll_seconds <= 1,
        "Invalid finalization wait",
    )
    started, session, stopped = time.monotonic(), None, False
    while True:
        current = status(native_root)
        if current.get("session_id") is not None:
            require(session in (None, current["session_id"]), "Worker session changed")
            session = current["session_id"]
        if current["state"] in ("NOT_STARTED", "COMPLETED", "FAILED", "CANCELLED"):
            verify(native_root)
            return {
                "terminal_acknowledged": True,
                "stop_requested": stopped,
                "wait_seconds": time.monotonic() - started,
                "owned_child_cleanup_verified": current.get(
                    "owned_child_cleanup_verified", current["state"] == "NOT_STARTED"
                ),
                "worker": current,
            }
        if time.monotonic() - started >= timeout:
            if session is not None and not stopped:
                stopped = finish(native_root, reason=reason)["state"] == "STOP_REQUESTED"
            return {
                "terminal_acknowledged": False,
                "stop_requested": stopped,
                "wait_seconds": time.monotonic() - started,
                "owned_child_cleanup_verified": False,
                "worker": current,
            }
        if request_stop and session is not None and not stopped:
            stopped = finish(native_root, reason=reason)["state"] == "STOP_REQUESTED"
        time.sleep(poll_seconds)


LAUNCHER = r"""@echo off
setlocal EnableExtensions DisableDelayedExpansion
echo Earth-Mercury: ENC search, local optimisation and native replay.
echo Leave this window open to follow the search and local optimisation.
echo RemoteSigned applies to this PowerShell process only. No permanent Windows setting changes.
set "WORKFLOW_SCRIPT=%~dp0RUN_NATIVE_WORKER.ps1"
powershell.exe -NoLogo -NoProfile -NonInteractive -Command "if ((Get-FileHash -LiteralPath $env:WORKFLOW_SCRIPT -Algorithm SHA256).Hash.ToLowerInvariant() -cne '@WORKER_SHA@') { throw 'Worker changed' }; $tokens=$null; $errors=$null; [System.Management.Automation.Language.Parser]::ParseFile($env:WORKFLOW_SCRIPT,[ref]$tokens,[ref]$errors) | Out-Null; if ($errors.Count -gt 0) { $errors | ForEach-Object { [Console]::Error.WriteLine($_.Message) }; exit 1 }; exit 0"
if errorlevel 1 goto failed
powershell.exe -NoLogo -NoProfile -NonInteractive -ExecutionPolicy RemoteSigned -File "%WORKFLOW_SCRIPT%"
set "WORKFLOW_EXIT=%ERRORLEVEL%"
echo Fixed worker ended. Exit code: %WORKFLOW_EXIT%
pause
exit /b %WORKFLOW_EXIT%
:failed
echo STOP: Worker identity or syntax check failed. No native process started.
pause
exit /b 1
"""
