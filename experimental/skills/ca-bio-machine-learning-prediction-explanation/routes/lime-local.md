# LIME for one prediction

```bash
python scripts/lime_local.py data.csv --label label --row 0 --positive AML --seeds 5 --out lime_results.tsv
```

Fits a random forest on a 70% split (refusing a test AUC under 0.6), then explains row `--row` with LIME under `--seeds` seeds (at least 2) and writes the local weights per seed plus the number of distinct top-k sets and the mean pairwise overlap. Options: `--num-features` (5), `--no-discretize`, `--kernel-width`.

- Report the overlap with the explanation. Rankings can change with seed, kernel width and discretization, especially in the tail; the head may be stable on a particular model. Measure it, do not assume either way.
- LIME fits a sparse linear surrogate on perturbed samples that fall off the data manifold for correlated genes. It is local: never average it into a global ranking.
- A biased model can be wrapped to look innocuous on the out-of-distribution points LIME and KernelSHAP probe, so a clean LIME explanation is not evidence a model is unbiased.
- For a local explanation of a tree model prefer the SHAP waterfall, labelled model-internal and background-relative.

Done when the explanation is reported with its top-k overlap across seeds and labelled local and model-internal.
