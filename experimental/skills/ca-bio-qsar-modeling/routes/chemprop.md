# Chemprop training

```bash
bash scripts/chemprop_pipeline.sh data.csv model_dir classification
```

Arguments are data path, output directory and classification or regression. The source defaults to five scaffold-split replicates, five models each, 50 epochs and seed 42; confirm installed 2.x flags before use. For prediction checkpoint paths and ensemble variance see `references/chemprop.md`. Split selection and applicability caveats are in `references/validation.md` and `references/selection.md`.

Done when models, split assignments and held-out performance are recorded.
