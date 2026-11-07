# FINEMAP or DAP-G

```bash
finemap --sss --in-files locus.master --n-causal-snps 5 --prob-conv-sss-tol 0.001 --n-iter 100000 --n-conv-sss 1000
```

Build the .z, signed space-delimited .ld and master files in identical SNP order using `references/finemap-cli.md`; `scripts/finemap_pipeline.sh` is an editable configuration example. For DAP-G use `references/dap-g-cli.md` and `scripts/dapg_finemap.sh`; the source build may return 1 even on success, so inspect output. Set the locus and ancestry instead of the example defaults.

Done when independent results and their alignment and convergence checks are reported.
