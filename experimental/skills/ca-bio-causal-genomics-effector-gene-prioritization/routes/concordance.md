# Integrate evidence

```bash
Rscript scripts/concordance_scoring.R locus_candidates.tsv out_prefix
```

Provide one row per locus and gene with the columns specified in the script header. Missing values do not pass and availability is reported separately. The source scores fine-mapping PIP > 0.5 and purity > 0.5, coloc >= 0.7, TSS distance <= 100 kb, PoPS top decile, L2G >= 0.5, and ABC >= 0.02 or rE2G >= 0.5; >=3 of six gives the source high tier. These thresholds are operational conventions, not proof of causality. Inspect overlap among evidence sources and method-specific calibrated thresholds; do not assume all streams independent. Use `references/evidence-framework.md` and `references/failure-modes.md`.

Done when scores, per-stream availability, discordance and uncertainty are reported.
