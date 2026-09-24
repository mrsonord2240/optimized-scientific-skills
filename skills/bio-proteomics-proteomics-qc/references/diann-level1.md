# Level-1 run metrics from a DIA-NN report

Read when a DIA-NN report is the only instrument-level evidence (RT fit, peak width).

## Level-1 Run Metrics From a DIA-NN Report

**Goal:** Read retention-time fit and peak width per run from the columns a DIA-NN report already carries.

**Approach:** Filter to Global.Q.Value <= 0.01, then per run correlate `RT` with `Predicted.RT` (R^2) and take the median `FWHM` (minutes in DIA-NN 2.x) and `Quantity.Quality`. With no rolling baseline yet, the across-run median stands in for it: a run whose FWHM is >25% above that median (the Thresholds row's 20-30% alarm) or whose RT fit R^2 is below 0.99 is flagged.

```python
import sys; sys.path.insert(0, "scripts")
from diann_level1 import diann_level1  # scripts/diann_level1.py
out = diann_level1(report)  # or: python scripts/diann_level1.py report.parquet
```

Three runs are a weak baseline: one broad run among three shifts the median. Trend these numbers over the queue (Levey-Jennings) before acting on a single flag.
