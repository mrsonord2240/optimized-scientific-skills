## Specificity Test (CellPhoneDB v5)

**Goal:** Get permutation specificity p-values with rigorous multi-subunit complex handling (human).

**Approach:** Run the statistical method on log-normalized counts plus a cell-type meta table; the permutation null shuffles cluster labels, and complexes require all subunits via the limiting (minimum) subunit.

```bash
python scripts/cellphonedb_statistical.py cellphonedb.zip meta.tsv counts_normalized.h5ad \
    --out-dir cpdb_out --iterations 1000 --threads 1 --seed 1337
```

The script's comments carry the parameter rationale, and its `__main__` guard is required on Windows: `score_interactions=True` uses `multiprocessing.Pool` internally, and without the guard the call crashes with `RuntimeError`. `debug_seed` fixes the permutation RNG (default `-1` is unseeded), but only at `threads=1`: with `threads>1` the pool workers do not replay the RNG stream (two identical seeded runs at `threads=4` differed on 82/120,375 significance flags). For a DEG-driven escape from one-vs-rest, use `cpdb_degs_analysis_method.call(..., degs_file_path=...)`.
