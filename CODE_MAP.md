# Where I start reading

Python files are in `workflow/code/local_optimizer/`. This is the fixed
Earth–Mercury ideal-sail case.

| Start here | What it does |
|---|---|
| `model.py` | Sail equations, RK4 segments and the shooting joins. |
| `optimiser.py` | Restoration or minimum-time problem, then the Ipopt call. |
| `validation.py` | Check the returned schedule, including one continuous flight from departure. |

For 600 holds: 3,600 later-state values + 1,200 angles + one time = **4,801 entries**.
Departure stays fixed. Restoration also fixes T at 600 days, so 4,800 entries
are free. Time refinement makes T free too.

`build_segment()` integrates the ideal-sail equations with RK4, using 16 substeps
per hold. CasADi builds the equations and their derivatives. Ipopt uses them to
adjust the angles, states and flight time. One Ipopt iteration is an update of
the proposed flight; the RK4 model performs the integration.

`model.build_problem()` joins the segments. `optimiser.build_problem()` adds the
phase objective and limits. For example, if a segment reaches 99 km and the next
guess says 100 km, the continuity constraint requires that 1 km gap to close.
The final commands are then propagated continuously from departure.

The 601 boundary states describe 600 command holds. Each hold lasts T/600 days.
Changing T changes the timing of all holds and Mercury's arrival state.

The connecting files:

- `ephemeris.py`: Mercury at the new arrival time.
- `refinement.py`: keep a checked result, enter an earlier time cell, or stop.
- `objective.py` and `controls.py`: arrival/time requirements and command-change limit.
- `case.py` and `workflow.py`: setup and calculation order.
- `bridge.py`, `native_results.py`, `transfer.py`: native jobs and exchanged data.
- `propagation.py`: independent continuous RK4 calculation.
- `records.py`, `console.py`, `report.py`: files, progress and saved-data plots.

`workflow/run.py` is the entry. `native_worker.ps1` handles the Windows jobs.
For native changes, read `workflow/InTrance_modified/InTrance/src/CInTrance.cpp`
and `CSimulation.cpp`. The bridge returns physical commands, not trained NN weights.

[Caruso](https://eprints.gla.ac.uk/226178/) gives the solar-sail multiple-shooting
example; their guess is shape-based, ours comes from ENC.
[CasADi](https://doi.org/10.1007/s12532-018-0139-4) supplies the derivatives;
[Ipopt](https://optimization-online.org/2004/03/836/) solves the constrained problem.
