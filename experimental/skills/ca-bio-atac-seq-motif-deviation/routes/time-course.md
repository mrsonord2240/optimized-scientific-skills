# Bulk time course

Run `routes/bulk.md` without `samples=` to get `chromvar_deviations.csv` (motifs x samples z-scores), then fit a spline regression of each motif's z-scores on time (5 or more time points) and look for non-monotonic trajectories. This step is not scripted and was not executed in testing.

Done when the answer lists motifs with significant time effects and the model used.
