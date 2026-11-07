# MultiVI (mosaic RNA + ATAC)

```bash
python scripts/multivi_mosaic.py mosaic.h5mu multivi_out
```

`mosaic.h5mu` holds `mod['rna']` (all cells, real counts) and `mod['atac']` (real counts for paired cells, all-zero rows for RNA-only cells, and the reverse for an ATAC-only block), both over the full cell union. MultiVI detects presence per cell from whether a block's counts sum to zero. Output: `multivi_out.h5mu` with `obsm['X_multivi']`. `--max-epochs N` overrides the default.

- Register with `setup_mudata`, not `setup_anndata` (deprecated for MULTIVI since scvi-tools 1.4; it warns, skips registration and `MULTIVI(adata)` then raises).
- No `batch_key`: the zero-fill mask is not a sequencing batch. Add `batch_key` only for real sequencing batches, as a separate covariate.
- Keep `--max-epochs` at the default for real runs: on 195 synthetic cells 40 epochs left RNA-only cells at chance placement (29% nearest paired same type) against 91% at the default.
- The missing modality is imputed, not measured: flag DE on it as model-dependent.
- The seed is fixed (`scvi.settings.seed = 0`).

Done when `multivi_out.h5mu` exists.
