# Report card and MultiQC

```bash
python scripts/aggregate_qc.py metrics.json --output atac_qc_mqc.tsv --sample NAME
```

`metrics.json` is `{metric: value}` for one sample (named by `--sample`) or `{sample: {metric: value}}`. Metric keys: `nuclear_reads_M`, `mt_fraction`, `NRF`, `PBC1`, `PBC2`, `TSS_enrichment`, `FRiP`. Missing, null or NaN values grade `NA` and the overall grade is `INCOMPLETE`; extra keys are ignored with a stderr note; a non-numeric value stops the run (exit 2). The output is a MultiQC custom-content TSV (keep the `_mqc.tsv` name) with the seven values, a grade per metric, and `overall_grade` (FAIL, then INCOMPLETE, then WARN, then PASS).

MultiQC ingests Picard insert-size metrics, samtools flagstat, deepTools fingerprint output and MACS peak tables. It does not compute TSS enrichment or NRF, so pass those through the `_mqc.tsv`:

```bash
multiqc fastqc/ picard/ samtools_stats/ macs2/ deeptools/ atac_qc_mqc.tsv -o multiqc_report
```

- Fill the JSON from `routes/complexity.md`, `routes/tss-encode.md`, `routes/alignment-mt.md` and `routes/fragment-tss-frip.md`.
- Drop or keep each replicate from the grades; a FAIL on TSS or FRiP is a library problem (`routes/diagnose.md`), not a reason to loosen the threshold.

Done when `atac_qc_mqc.tsv` exists and each sample has a keep or drop decision.
