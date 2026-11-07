# Paired observations

```bash
Rscript scripts/annotate_paired.R paired.csv paired.png holm
```

CSV requires `subject_id,time,value`, exactly two levels and one complete row per subject and level. The script verifies matching IDs before the paired Wilcoxon test and writes an image and results CSV. Python alternative `python scripts/annotate_paired.py paired.csv paired.png holm` performs the same ID check; use textual time levels because its pivot uses original column types after deriving string labels. The single comparison's correction does not account for comparisons outside this invocation. A paired effect magnitude and uncertainty are not supplied; calculate and report a design-appropriate measure separately. Outputs may be overwritten.

Done when subject matching, contrast direction, exact p-value, family and effect magnitude are documented.
