# scCODA composition
~~~sh
python scripts/sccoda_run.py annotated.h5ad credible_effects.csv
~~~
Input obs needs sample, condition and cell_type; each sample must belong to one condition with independent replicates. Source settings use automatic reference, TensorFlow seed 42, default HMC and estimated FDR 0.1. All effects are relative to the selected reference; assess its stability, convergence and sensitivity to alternative plausible references. Verify installed API and dependencies.
**Done when:** credible effects, reference, diagnostics, settings and sample counts are reviewed.
