# Many traits

```bash
Rscript scripts/hyprcoloc_traits.R betas.tsv ses.tsv hyprcoloc_results.tsv
```

Each input is a tab-separated SNP-by-trait matrix with SNP IDs in the first column and matching trait headers; rows and columns must be identical and ordered alike. The wrapper uses the source thresholds reg.thresh=0.7 and align.thresh=0.7. HyPrColoc assumes one causal variant per cluster; see `references/hyprcoloc.md`. For 3-5 traits, `references/moloc.md` supplies an input-dependent example; no file-based moloc wrapper is provided, so adapt it to validated trait tables and the installed interface.

Done when trait clusters or multi-omic posteriors and assumptions are reported.
