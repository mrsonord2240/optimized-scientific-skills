# scVI and scANVI details

Use [`../scripts/scvi_scanvi_integration.py`](../scripts/scvi_scanvi_integration.py) when the input is an `.h5ad` with raw counts in `X`. The script copies them to `layers['counts']` before HVG selection, sets `scvi.settings.seed`, and takes an explicit `--max-epochs` value so a run can be reproduced.

For scANVI, the label column must contain the literal `Unknown` for unlabeled cells. `SCANVI.from_scvi_model(model, 'Unknown', labels_key=...)` takes the unlabeled category as its second positional argument. Interpret `X_scVI` and `X_scANVI` only as embeddings for neighbours, clusters, and visualization.
