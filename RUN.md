# Running the workflow

Windows runs the existing InTrance build. I developed the local solver and plots
in Python on Mac, so I connected the two through a shared folder. Both parts use
the same case files and pass the results back for checking.

## Set up once

1. Open this repository's **Releases** page. Download the attached
   `intrance-local-workflow.zip` and extract it into a folder that Mac and
   Windows can both see. This ZIP includes the compiled InTrance EXE.
   Use the attached ZIP; GitHub's source-code archive does not include the EXE.
2. Put the four supplied data files in `native_data/`:
   `bin1600_600.405`, `ELEMENTS.UNNUM`, `ELEMENTS.NUMBR`, `ELEMENTS.COMET`.
3. In Mac Terminal, change to the extracted folder. With Python 3.12 installed:

```sh
python3.12 -m venv .venv-mac
.venv-mac/bin/python -m pip install -r workflow/requirements.txt
```

Create this environment after putting the folder at its final location.

4. Open this same folder in Windows. Run `CHECK_NATIVE.cmd` once.
   It should finish with `Native interface checked.` and exit code 0.
5. Back on Mac:

```sh
sh START_MAC.command check --run-name demo_01
```

Expected: `CHECKED_NO_SCIENCE`. Setup is ready.

## Run

On Mac:

```sh
caffeinate -i sh START_MAC.command run --run-name demo_01
```

When Mac says `Waiting for the Windows worker`, open
`START_WINDOWS_WORKER.cmd` in the same folder on Windows. Keep both terminals
open. The files move between stages automatically.

```text
Windows: ENC search -> export -> replay starting commands
Mac:     fixed-time restoration -> time refinement -> checks
Windows: replay the selected schedule from departure
Mac:     final comparison -> report and plots
```

The default uses seed 1 and 200 completed ENC evaluations to obtain the starting
controller. Its full 600-day history provides 600 held commands and 601 boundary
states. Evaluations count search work; holds divide one flight.

Each local solve gets one console header. `TOF now [d]` is the current time being tried;
`Arrival cost` is used while time stays fixed. Constraint and dual errors help
read the solve. The returned schedule then goes through the flight checks.

Wait for `Workflow passed.` and then `Run complete.` in the Mac window.
A complete run in the tested setup took about 13 min 31 s, including waits and
checks. Runtime depends on the machine and worker readiness. Finish setup
before the demonstration and explain the workflow while the calculation runs.

Run one attempt at a time in this folder. The Windows launcher follows the run
selected by Mac. Start it once, after the waiting message.

Use a new name for another attempt, e.g. `demo_02`. Keep a stopped run's files
so the reason can be checked.

## Limits used here

| Setting | Value |
|---|---|
| Local phase allowance | 3,000 s, including local solves and checks |
| One Ipopt call | Up to 900 s and 3,000 iterations; less time if the budget is nearly used |
| Attempts in one time cell | Up to 3 in total; retries need a checked, improving result stopped by a resource limit |
| Native search / one export or replay | 720 s / 90 s |
| Windows worker session | 7,200 s, including waiting for Mac |

A new local call needs at least 180 s left. Its process allowance is at most
1,200 s, clipped to the remaining budget. These are practical settings I chose.

The recorded result used five time solves and four handovers between arrival-time
cells. A checked result at a lower cell boundary can enter the next earlier cell
while the budget permits. Convergence inside a cell ends that start. The declared
flight-time range is 0.5--600 days.

## Read the output

Open `runs/shared/demo_01/report/README.md`. Main numbers are also in
`result_summary.json`.

| Figure | What I check |
|---|---|
| `trajectory_quality` | Native paths and the two reporting scores |
| `controls` | Clock/cone commands and consecutive normal changes |
| `refinement_progress` | Kept flight times and the last Ipopt log |
| `numerical_checks` | Native, finer RK4 and shooting-node differences |

Each figure has PDF, PNG and SVG. `plot_index.json` lists what was generated.
RK4-128 and RK4-256 mean 128 and 256 substeps per command hold. The native
comparison uses 100 m / 0.001 m/s limits; the finer-grid check uses 1 m /
0.00001 m/s. Arrival errors are checked separately against the reference caps.

## Show a saved result

For the result included in this download:

```sh
.venv-mac/bin/python workflow/run.py plot-saved \
  --result results/mercury_561d --output-dir plots/saved_561d
```

For all four views from a complete run folder:

```sh
.venv-mac/bin/python workflow/run.py plot runs/shared/demo_01 \
  --output-dir plots/demo_01
```

These commands read saved data. Choose a new output folder when redrawing.
The short `results/mercury_561d` selection contains the trajectory and commands;
full solver logs and continuous-check arrays stay with the complete run.
