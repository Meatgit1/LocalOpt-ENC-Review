"""Read the evaluated controller and its ENC trace before using it as a local start."""

import csv
import math
import re
from pathlib import Path

from .records import require, sha, load, input_file


BANK_OUTPUTS = ("final.eac", "final.json", "summary.json", "trace.csv", "timing.csv")


def search_contract(
    seed, evaluation_limit, search_seconds, terminal_policy="exact_completed_count"
):
    require(
        type(seed) is int and 0 <= seed <= 2147483647,
        "Explicit native seed must be a nonnegative integer",
    )
    require(
        type(evaluation_limit) is int
        and 40 <= evaluation_limit <= 1000000
        and evaluation_limit % 2 == 0,
        "Cutoff must be an even returned-evaluation count between 40 and 1000000",
    )
    require(
        type(search_seconds) is int and 1 <= search_seconds <= 1200,
        "Search child deadline must be 1..1200 seconds",
    )
    require(
        terminal_policy == "exact_completed_count",
        "Review workflow requires an exact completed-evaluation cutoff",
    )
    return {
        "seed": seed,
        "evaluation_limit": evaluation_limit,
        "child_seconds": search_seconds,
        "stop_policy": terminal_policy,
        "warm_start": False,
        "deadline_policy": "hard_timeout_reject_not_orderly_cutoff",
    }


def _eac(path):
    tokens = Path(path).read_text().split()
    require(len(tokens) == 882 and tokens[0] == "881", "Wrong controller dimension")
    values = [float(x) for x in tokens[1:]]
    require(all(math.isfinite(x) for x in values), "Nonfinite controller")
    return values


def _fnv(data):
    value = 2166136261
    for b in data:
        value = ((value ^ b) * 16777619) & 0xFFFFFFFF
    return f"{value:08x}"


def inspect_observer(result, search):
    """Check the complete evaluation count and the controller that actually produced it."""
    result = Path(result)
    require(
        search
        == search_contract(
            search["seed"],
            search["evaluation_limit"],
            search["child_seconds"],
            search["stop_policy"],
        ),
        "Invalid search contract",
    )
    bank = result / "enc_observer"
    require({p.name for p in bank.iterdir()} == set(BANK_OUTPUTS), "Unexpected observer bank files")
    summary, final = load(input_file(bank, "summary.json")), load(input_file(bank, "final.json"))
    require(
        set(summary)
        == {
            "schema",
            "status",
            "study",
            "seed",
            "capture",
            "evaluation_limit",
            "completed_evaluations",
            "completed_populations",
            "completed_pairs",
            "stop_reason",
            "controller_fnv1a32",
            "checkpoints",
            "timing_semantics",
        }
        and type(summary["schema"]) is int
        and summary["schema"] == 1
        and summary["status"] == "complete"
        and summary["study"] == "count_bank_only"
        and summary["capture"] is False
        and summary["seed"] == search["seed"]
        and summary["evaluation_limit"] == search["evaluation_limit"]
        and summary["timing_semantics"]
        == "process_clock_ticks_platform_dependent_not_campaign_wall_time",
        "Observer cutoff/seed/stop contract differs",
    )
    require(
        all(
            type(summary[k]) is int
            for k in (
                "seed",
                "evaluation_limit",
                "completed_evaluations",
                "completed_populations",
                "completed_pairs",
            )
        ),
        "Observer counters must be integers",
    )
    count, cap, reason = (
        summary["completed_evaluations"],
        search["evaluation_limit"],
        summary["stop_reason"],
    )
    require(40 <= count <= cap and count % 2 == 0, "Invalid attained count")
    require(
        reason == "evaluation_limit" and count == cap,
        "Observer cutoff/seed/stop contract differs",
    )
    actual = cap if count == cap else None
    require(
        summary["checkpoints"] == [{"requested": cap, "actual": actual, "captured": False}]
        and type(summary["checkpoints"][0]["requested"]) is int
        and (
            summary["checkpoints"][0]["actual"] is None
            if actual is None
            else type(summary["checkpoints"][0]["actual"]) is int
        )
        and summary["checkpoints"][0]["captured"] is False,
        "Observer checkpoint types differ",
    )
    with input_file(bank, "trace.csv").open(newline="") as stream:
        reader = csv.DictReader(stream)
        require(
            reader.fieldnames
            == [
                "completed_evaluations",
                "completed_populations",
                "completed_pairs",
                "scan_epoch",
                "refinement_epoch",
                "hypercube",
                "mutation_scale",
                "best_fitness",
                "controller_fnv1a32",
                "event",
            ],
            "Wrong observer trace columns",
        )
        rows = list(reader)
    require(rows and rows[0]["event"] == "population", "Observer initial population absent")
    count = populations = pairs = 0
    for row in rows:
        require(
            None not in row
            and all(v is not None for v in row.values())
            and row["event"] in ("population", "pair"),
            "Partial observer trace",
        )
        population = row["event"] == "population"
        count += 40 if population else 2
        populations += int(population)
        pairs += int(not population)
        require(
            (
                int(row["completed_evaluations"]),
                int(row["completed_populations"]),
                int(row["completed_pairs"]),
            )
            == (count, populations, pairs)
            and count <= search["evaluation_limit"]
            and re.fullmatch(r"[a-f0-9]{8}", row["controller_fnv1a32"]),
            "Observer count trace differs",
        )
        require(
            all(
                math.isfinite(float(row[k]))
                for k in ("hypercube", "mutation_scale", "best_fitness")
            )
            and all(int(row[k]) >= 0 for k in ("scan_epoch", "refinement_epoch")),
            "Invalid observer trace metrics",
        )
    require(
        (count, populations, pairs)
        == (
            summary["completed_evaluations"],
            summary["completed_populations"],
            summary["completed_pairs"],
        ),
        "Observer summary differs from trace",
    )
    expected_keys = {
        "schema",
        "seed",
        "reason",
        "requested_evaluations",
        "completed_evaluations",
        "completed_populations",
        "completed_pairs",
        "scan_epoch",
        "refinement_epoch",
        "hypercube",
        "mutation_scale",
        "best_fitness",
        "parameter_count",
        "controller_fnv1a32",
    }
    require(
        set(final) == expected_keys
        and type(final["schema"]) is int
        and final["schema"] == 1
        and type(final["seed"]) is int
        and final["seed"] == search["seed"]
        and final["reason"] == reason
        and type(final["parameter_count"]) is int
        and final["parameter_count"] == 881
        and type(final["requested_evaluations"]) is int
        and final["requested_evaluations"] == final["completed_evaluations"] == count,
        "Observer final metadata differs",
    )
    for key in (
        "completed_evaluations",
        "completed_populations",
        "completed_pairs",
        "scan_epoch",
        "refinement_epoch",
    ):
        require(
            type(final[key]) is int and final[key] == int(rows[-1][key]),
            "Final controller count/epoch differs",
        )
    for key in ("hypercube", "mutation_scale", "best_fitness"):
        require(
            type(final[key]) in (int, float) and math.isfinite(final[key]),
            "Nonfinite final observer metric",
        )
    # The final controller must match the last evaluated trace record.
    metric_keys = ("hypercube", "mutation_scale", "best_fitness")
    for key in metric_keys:
        require(
            type(final[key]) in (int, float)
            and math.isfinite(final[key])
            and final[key] == float(rows[-1][key]),
            "Final observer metrics differ",
        )
    digest = _fnv(input_file(bank, "final.eac").read_bytes())
    require(
        digest
        == final["controller_fnv1a32"]
        == summary["controller_fnv1a32"]
        == rows[-1]["controller_fnv1a32"],
        "Observer final controller differs from trace",
    )
    values = _eac(bank / "final.eac")
    require(
        _eac(input_file(result, "mercury.eac")) == [float(format(x, ".16g")) for x in values],
        "Saved16/observer17 controller differs",
    )
    with input_file(bank, "timing.csv").open(newline="") as stream:
        reader = csv.DictReader(stream)
        require(
            reader.fieldnames
            == ["completed_evaluations", "process_clock_ticks", "clocks_per_second"],
            "Wrong observer timing columns",
        )
        timing = list(reader)
    require(len(timing) == len(rows), "Incomplete observer timing")
    previous = -1
    for row, trace in zip(timing, rows):
        ticks = int(row["process_clock_ticks"])
        require(
            row["completed_evaluations"] == trace["completed_evaluations"]
            and ticks >= previous
            and int(row["clocks_per_second"]) > 0,
            "Invalid observer timing",
        )
        previous = ticks
    inspected = {
        "seed": search["seed"],
        "completed_evaluations": count,
        "completed_populations": populations,
        "completed_pairs": pairs,
        "stop_reason": reason,
        "best_fitness": final["best_fitness"],
        "controller_sha256": sha(bank / "final.eac"),
        "native_selected_endpoint": (
            "retained in mercury.csv/mercury.rep; not the full export " "horizon"
        ),
    }
    return inspected


def selected_endpoint(result):
    """Read the native selected endpoint's printed physical values, not subfitness.

    Native CSV is rounded output. These numbers are descriptive search metrics;
    they are not a high-precision independent replay or final acceptance gate.
    """
    path = input_file(Path(result), "mercury.csv")
    columns = ("TOF [d]", "MJD [d]", "Dist (MERCURY) [km]", "Vel (MERCURY) [km/s]")
    with path.open(newline="") as stream:
        reader = csv.DictReader(stream, skipinitialspace=True)
        require(
            reader.fieldnames is not None and all(reader.fieldnames.count(n) == 1 for n in columns),
            "Native selected endpoint CSV columns absent/repeated",
        )
        previous, count, last = -1.0, 0, None
        for row in reader:
            require(
                None not in row and all(row.get(n) is not None for n in columns),
                "Partial native CSV row",
            )
            values = [float(row[n]) for n in columns]
            require(
                all(math.isfinite(v) for v in values)
                and values[0] > previous
                and values[2] >= 0
                and values[3] >= 0,
                "Invalid native selected endpoint history",
            )
            if count == 0:
                require(values[0] == 0, "Native CSV must begin at departure")
            previous, last, count = values[0], values, count + 1
    require(count >= 2, "Native selected endpoint history is empty")
    return {
        "flight_time_days": last[0],
        "arrival_mjd_days": last[1],
        "position_miss_km": last[2],
        "velocity_miss_m_s": last[3] * 1000.0,
        "rows": count,
        "basis": (
            "last row of native-selected mercury.csv; printed precision, " "not full-horizon export"
        ),
        "high_precision_acceptance": False,
    }
