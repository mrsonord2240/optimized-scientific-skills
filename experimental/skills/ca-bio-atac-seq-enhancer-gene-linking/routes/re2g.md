# ENCODE-rE2G

```bash
snakemake -j1 --use-conda
```

Run from the configured ENCODE-rE2G checkout. Set config/config.yaml and the biosamples TSV with accessibility, optional H3K27ac and matched or megamap Hi-C. The pipeline selects a model from input types; it has no no-Hi-C/powerlaw model. Use its model-specific threshold and disclose proxy contacts. See `references/method-reference.md` and `references/setup.md`.

Done when calibrated predictions and model and contact provenance are reported.
