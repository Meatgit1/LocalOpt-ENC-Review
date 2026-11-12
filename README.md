# InTrance + local optimisation

This repository contains the Earth–Mercury workflow, selected source and
recorded results. Full project:
[enc-rendezvous-thesis](https://github.com/Meatgit1/enc-rendezvous-thesis).

Start on Mac. Windows runs InTrance to get a starting flight. Mac changes the
sail commands and time. Windows replays that schedule from departure; Mac
checks the result and makes the plots.

The recorded flight is **561.885720 days**, compared with the **576-day** native
reference. Arrival misses: **62,989.754 km / 63.024693 m/s**, inside the
reference-quality caps.

[Download the runnable ZIP](https://github.com/Meatgit1/LocalOpt-ENC-Review/releases/latest/download/intrance-local-workflow.zip)

The ZIP includes the native EXE. Put it in a shared folder and follow
[RUN.md](RUN.md). The separate data files and Python environment are needed once.

- [RUN.md](RUN.md): setup and commands.
- [CODE_MAP.md](CODE_MAP.md): where to start reading.
- [561.886-day result](results/mercury_561d/README.md): settings, checks and saved files.
- [Parallel test](results/parallel_search/README.md): what I tried and what is still open.

Entry: `workflow/run.py`. Python code: `workflow/code/local_optimizer/`.
Native source: `workflow/InTrance_modified/`.

Copied source: [main project snapshot](https://github.com/Meatgit1/enc-rendezvous-thesis/commit/a712b84298a5ec9965e58a8006006b7e895d775f).

InTrance is existing work by Bernd Dachwald and Andreas Ohndorf. Its original
notices remain in the source and [NOTICE.md](NOTICE.md).
