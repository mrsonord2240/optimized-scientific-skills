# Peptide or precursor table

Two-group test from a MaxQuant `evidence.txt` with msqrob2:

```bash
Rscript scripts/msqrob2_de.R evidence.txt annotation.csv Treatment-Control results.csv
```

`annotation.csv`: `Raw.file, Condition, BioReplicate, IsotopeLabelType`. The script writes `results.csv` (`logFC, se, df, t, pval, adjPval`) and `results_undetected.csv`, and stops when the median logFC exceeds 0.05.

- `model=peptide` keeps every peptide as its own degree of freedom (`~condition + (1|sample) + (1|feature)`). Use it when peptides within a protein disagree or coverage is unbalanced.
- msqrob2 `ridge = TRUE` is refused on a two-group design. Use it only for three or more groups or group plus covariates.
- Undetected proteins (no observation in a condition) have `adjPval = NA`. They are in `results_undetected.csv`.

Technical replicates, nested or repeated measures, SRM, PRM or DIA: MSstats.

```bash
Rscript scripts/msstats_group_comparison.R evidence.txt proteinGroups.txt annotation.csv Control Treatment FALSE out_prefix
```

- Writes `out_prefix_tested.csv` (`adj.pvalue` is the BH p) and `out_prefix_undetected.csv` (`issue == oneConditionMissing`, infinite `log2FC`).
- The sixth argument is the normalization. `FALSE` normalizes nothing. `equalizeMedians` (the MSstats default) centres each run on its own detected peptides and makes the test anticonservative unless the centring checks pass.
- Normalize on features present in every run, or at protein level after summarization.

Centring checks when you normalize the peptide table: `references/feature_level.md`; the residual-SD and offset checks are in `scripts/centring_checks.R`.

Done when the tested table and the undetected list exist and the median logFC is near 0.
