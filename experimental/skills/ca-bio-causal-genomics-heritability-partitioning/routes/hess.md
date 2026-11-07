# Local heritability

```bash
bash scripts/hess_local_h2.sh sumstats.txt 1 ref_chr1 partition_chr1.bed hess
```

Set HESS_DIR to a compatible source checkout; the source Python 3 patch requirements and input details are in `references/hess-local-h2.md`. Repeat step 1 for all autosomes before the separate aggregate step; that step was not exercised by the source. Flag sparse loci.

Done when per-locus estimates and uncertainty are reported with reference and coverage.
