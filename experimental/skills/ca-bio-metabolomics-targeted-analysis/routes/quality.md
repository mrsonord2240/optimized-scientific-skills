# Assay quality summaries

```bash
Rscript scripts/precision_nested_anova.R file=qc.csv
```

Input CSV needs `qc_level,day,measured_conc`; the script computes within-day and balanced one-way variance-component precision and requires at least two days with equal replicate counts. `level=LLOQ` restricts the level. Other source checks are `Rscript scripts/matrix_recovery_carryover.R matrix_factor file=lots.csv` and `Rscript scripts/matrix_recovery_carryover.R recovery file=recovery.csv`; inspect their headers for required columns and arguments. Carryover mode expects a blank response ratio and a compatible ratio-calibration slope and intercept, and does not assess IS carryover unless separate inputs are supplied. These arithmetic outputs and example thresholds are not a complete assay validation; inspect units, denominators, replicate structure and the governing protocol. Preserve failure reasons rather than converting missing results to zero.

Done when precision, matrix, recovery, carryover and other required validation evidence have been reviewed in context.
