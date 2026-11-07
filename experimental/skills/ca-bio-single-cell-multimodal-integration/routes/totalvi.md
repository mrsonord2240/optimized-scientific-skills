# totalVI (CITE-seq)

```bash
python scripts/totalvi_cite_seq.py cite_seq.h5mu totalvi_out
```

`cite_seq.h5mu` holds `mod['rna']` (raw counts in layer `counts`) and `mod['prot']` (raw ADT counts in `.X`). Output: `totalvi_out.h5mu` (`obsm['X_totalVI']`), `totalvi_out_foreground.csv` (per-protein foreground probability, 1 minus the background weight), `totalvi_out_denoised_prot.csv`. `--max-epochs N` overrides scvi-tools' dataset-size default.

- Training is stochastic: two unseeded runs differed by up to 0.97 in latent values. The script sets `scvi.settings.seed = 0`.
- totalVI models each protein as a two-component mixture (background and foreground), so no separate DSB is needed.
- It overfits tiny datasets.

Done when `totalvi_out.h5mu` exists.
