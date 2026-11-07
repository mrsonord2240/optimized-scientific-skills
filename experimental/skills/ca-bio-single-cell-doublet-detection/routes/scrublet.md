# Scrublet for one capture
~~~sh
python scripts/scrublet_run.py capture.h5ad scored.h5ad --expected-rate 0.04 --seed 0
~~~
Input X must contain raw counts after basic QC from one physical capture. Set the rate from the full capture loading and platform guidance; the source's rough 0.008 per 1000 recovered cells is a heuristic, not a universal calibration. Script saves scores and flags without deleting cells. Inspect score distribution and threshold; homotypic doublets and high-RNA true cells can be misclassified.
**Done when:** scored output, rate, threshold and uncertain cells are reviewed before any filtering.
