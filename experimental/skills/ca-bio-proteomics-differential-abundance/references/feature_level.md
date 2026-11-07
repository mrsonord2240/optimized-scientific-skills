# Feature-level testing: centring checks and MSstats detail

Read when a peptide table is normalized before msqrob2 or MSstats.

## Centring checks

Per-run median normalization of a peptide table gives an offset between conditions and a smaller residual. Run both checks before reading any feature-level result. `scripts/centring_checks.R` defines them:

```bash
Rscript scripts/centring_checks.R peptide_wide.csv samples.csv TRUE Control Treatment
```

`peptide_wide.csv`: `feature`, `protein`, one linear intensity column per run. `samples.csv`: `run`, `condition`. `TRUE` applies per-run median centring first.

1. Residual SD. msqrob2: after `filterFeatures`, before `aggregateFeatures`, compare the residual SD of the peptide table before and after normalizing; a ratio below 0.85 warns. MSstats: fit with `normalization = FALSE` and with the default, and compare `median(tested$SE)`; same ratio.
2. Offset. The median `log2FC` of the result table (msqrob2 `res$logFC`; MSstats `tested$log2FC`) above 0.05 in absolute value stops the analysis.

- A failure is fixed upstream: normalize on features present in every run, or at protein level after summarization. Do not re-centre p-values.
- Passing both checks does not certify the normalization. A study that expects a global shift sets `OFFSET_MAX` deliberately.
- Read a drop in residual SD together with the offset: real loading variation also shrinks it.

## MSstats

- `MaxQtoMSstatsFormat` imports, `dataProcess` summarizes (TMP), `groupComparison` tests.
- The contrast matrix is columns-by-condition; its `dimnames` must equal `levels(proc$ProteinLevelData$GROUP)`. The script builds it from the two level names.
- Keep `MBimpute = FALSE`. It is the only place MSstats imputes (AFT censored imputation); `groupComparison` has no censoring argument.
- `adj.pvalue` is the BH p.
