## CITE-seq: totalVI (Python, denoise + DE in one model)

**Goal:** Jointly model RNA + protein with explicit protein background, yielding a denoised latent space and foreground probabilities.

**Approach:** Register a MuData object, train the conditional VAE, then read the latent representation and per-protein foreground probability.

```bash
python scripts/totalvi_cite_seq.py cite_seq.h5mu totalvi_out   # -> totalvi_out.h5mu (obsm['X_totalVI']), totalvi_out_foreground.csv, totalvi_out_denoised_prot.csv
```

`cite_seq.h5mu` holds `mod['rna']` (raw counts in layer `counts`) and `mod['prot']` (raw ADT counts in `.X`). scvi-tools VAE training is stochastic unless seeded: two unseeded runs on identical input differed by up to 0.97 (max abs latent diff), and the script sets `scvi.settings.seed = 0` so reruns are bit-identical.

## Mosaic: MultiVI (Python, RNA+ATAC partially observed)

**Goal:** Jointly embed a mosaic design -- some cells have both RNA and ATAC (paired), others only one modality -- imputing the missing side.

**Approach:** Build one MuData with an RNA AnnData and an ATAC AnnData that both cover the full cell union; cells missing a modality get all-zero rows for that modality's block (MultiVI detects presence per cell from whether that block's raw counts sum to zero, not from a separate flag). Register with `setup_mudata`, not `setup_anndata` -- `MULTIVI.setup_anndata` on a plain AnnData is deprecated since scvi-tools 1.4 and silently skips registration (warns, then `MULTIVI(adata)` raises "Please set up your AnnData with MULTIVI.setup_anndata first").

```bash
python scripts/multivi_mosaic.py mosaic.h5mu multivi_out     # -> multivi_out.h5mu, obsm['X_multivi']
```

`mosaic.h5mu` holds `mod['rna']` (all cells, real counts) and `mod['atac']` (real counts for paired cells, all-zero rows for RNA-only cells, and vice versa for an ATAC-only block). The script seeds `scvi.settings.seed = 0`. Leave `--max-epochs` at scvi-tools' default for real runs: on 195 synthetic cells, 40 epochs left RNA-only cells at chance placement (29% nearest-paired same type) while the default run reached 91%.

Verified on synthetic 150-cell mosaic data (90 paired, 60 RNA-only, 3 known cell types, scvi-tools 1.5.1): RNA-only cells land nearer their same-type paired counterparts (mean latent distance 0.27) than different-type ones (0.56), confirming the model actually uses the shared RNA signal to place unpaired cells rather than clustering by modality of origin. This example passes no `batch_key` because the modality mask above is not a sequencing batch; if cells also span real sequencing batches, add `batch_key` for that separately -- see Common Errors' MultiVI row for the pitfall of confusing the two.
