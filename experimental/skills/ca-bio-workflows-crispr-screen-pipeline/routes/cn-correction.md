# Copy-number correction (cancer cell lines)

```bash
Rscript scripts/cn_correction.R experiment.count.txt screen
```

Writes `screen_cleanr_corrected_counts.txt`. Pass that file, not `experiment.count.txt`, as the count table to every hit caller. Its counts are non-integer and guides below `min_reads` (30, in the first sample column) are dropped, so `scripts/qc.py` refuses it: run QC on the raw table first.

- The built-in annotation is `KY_Library_v1.0`; for another library pass `library=<annotation.tsv>` (columns CODE, GENES, EXONE, CHRM, STRAND, STARTpos, ENDpos).
- The baseline sample is the first sample column.
- Amplified loci look essential whatever the gene does (ERBB2, MYC, FGFR1): correct before hit calling. CRISPRi (Dolcetto) avoids the cut; its knockdown is less complete.
- CRISPRcleanR needs no CN profile, so it cannot check itself. To verify, compute Spearman between gene LFC and CN with a matched profile: absolute rho below 0.10 after correction (0.05 stricter).
- Matched CN for many lines with longitudinal data: `routes/chronos.md` instead.
- Not a cancer line: skip this route.

Done when the corrected count file exists and the hit callers read it.
