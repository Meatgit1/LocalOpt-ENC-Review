# Trying two ENC starts together

Should I spend longer in ENC, or give the local optimiser an earlier guess?
I tried seeds **11 and 29**, running two independent chains together.
Same physics, 600 controls, arrival caps and command-change limit.

| Handover choice | Result for both seeds |
|---|---|
| After 200 evaluations | Local refinement reached **585.5 days** and passed native checks. |
| Up to 10,000 evaluations | Native search hit its **720 s** deadline. No final handover, so no local result. |

The early chains stopped when too little local budget remained for another solve.

Observed hypercube size stayed **H = 1**. The parallel work here was two
independent seeded searches.

The deeper searches did not reach local refinement. The best number of ENC
evaluations before handover is still open. Their last complete trace counts
were 4,204 and 4,136.

Next: save complete intermediate controllers during a longer search and keep
time available for their local refinement.

The measured outcomes are in `results.csv`. The normal launcher runs one seed;
this folder keeps the parallel test result for discussion.
