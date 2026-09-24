# Harmony details

Use [`../scripts/harmony_integration.py`](../scripts/harmony_integration.py) for the complete seeded preprocessing and integration path.

`theta` is a diversity penalty, not an on/off switch. Start near 2 and change it only after paired batch-mixing and biological-conservation metrics show under-correction. Values around 0.5--5 can be visually indistinguishable when the batch effect is modest; `theta = 0` still performs correction. Preserve the uncorrected PCA and compare rare populations after every change.
