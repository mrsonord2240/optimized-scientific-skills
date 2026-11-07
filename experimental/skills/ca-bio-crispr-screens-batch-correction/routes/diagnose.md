# Diagnose
~~~bash
python scripts/batch_diagnostic.py counts.tsv metadata.tsv > diagnostic.txt
~~~
Use feature-indexed tab-separated counts and sample-indexed metadata containing batch and condition. Check all count samples are in metadata; the source script silently uses overlaps. Inspect PCA, ANOVA and confounding; an F statistic alone is insufficient.
Done when diagnostic output and sample allocation have been reviewed.
