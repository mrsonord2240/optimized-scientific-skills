# Analyze multiome accessibility
~~~bash
python scripts/multiome_differential.py multiome.h5mu --group "GENE_A KO" --control NTC --out peaks.tsv
~~~
Use existing matched RNA and ATAC modalities with shared cells; RNA must contain per-target mixscape_class. Replace group and control with actual labels. The source uses normalized log1p counts and Wilcoxon, not TF-IDF, and requires independent evaluation of statistical assumptions. Peak-to-gene linking requires annotation; see [multiome notes](../references/multiomic-perturb-seq.md).
Done when cell alignment, controls and peak results are reviewed.
