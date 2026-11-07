# TreeSHAP with module aggregation

```bash
python scripts/shap_tree.py data.csv --label label --positive AML --mode both --out shap_results.tsv
```

`data.csv`: one row per sample, a binary label column, numeric feature columns (`--drop` for ids and batch). The script fits a 300-tree random forest on a 70% split, refuses a test AUC under 0.6, then runs TreeSHAP with `feature_perturbation` set explicitly. It builds a co-expression module map (|Spearman| on training rows, average linkage, cut at `--cut` 0.7 in 1 - |rho|) and writes one row per module and per gene with mean |SHAP|. More than `--max-features` (2000) features keeps the highest-variance ones, label-free.

- `--mode interventional` (default) answers "what does the model use": features replaced from a training background (`--background-n`, 100-1000 rows). `tree_path_dependent` is conditional and answers "what is informative here"; do not call it model reliance. `both` prints the top-20 gene overlap between them.
- A gene absent from every tree split gets zero under either mode. Correlated genes the trees do use can swap credit between modes; the module total is the stable quantity.
- Rank the modules, not genes inside one. The module map is an analysis choice: report the correlation, linkage and cut, and check sensitivity to the cut.
- The background defines "absence of a feature". Compare two defensible backgrounds (for example tumour-heavy and healthy-tissue) and report the top-k overlap; a single global mean of a heterogeneous cohort is no real sample.
- The additivity check stays on. If it fails, the cause is usually raw margin against probability; fix the configuration. Attributions in log-odds differ from probability space: state the scale.
- To explain a model you already trained, call `positive_class_values` and `module_map` from the script on your own explainer output. SHAP plotting: new-style `shap.plots.*` take an `Explanation`, legacy `summary_plot` takes arrays; do not mix them.

Done when `shap_results.tsv` exists and the answer ranks modules, states the mode, background and scale, and says the ranking is model-internal.
