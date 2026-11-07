# Tissue enrichment

```bash
ldsc.py --h2-cts trait.sumstats.gz --ref-ld-chr baseline. --ref-ld-chr-cts tissues.ldcts --w-ld-chr weights. --out trait_cts
```

The source CBIIT version requires a patch for h2-cts; see `references/compatibility.md`. Use matched chromosome-split references and the manifest contract in `references/cell-type-prioritization.md`. Correct across tested tissues.

Done when ranked tissue results and multiplicity are reported.
