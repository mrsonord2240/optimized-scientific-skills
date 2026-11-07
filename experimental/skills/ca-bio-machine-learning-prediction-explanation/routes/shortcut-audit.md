# Using attributions to catch shortcut and batch learning

No script; this is the most trustworthy use of attribution. Run `routes/shap-tree.md` first, then look at what the top modules are.

- If top features are batch indicators, sequencing-depth-tracking housekeeping genes, or other technical correlates, the model is cheating. Confirm with the batch checks in `ca-bio-machine-learning-omics-classifiers` (its batch-check route: batch predictability and leave-one-batch-out).
- A high AUC with technical top features means the explanation exposed an artifact; report that, not a gene list.
- If the audit is clean, the ranking is still model-internal. It is not causal (predictive features mix direct effects, confounders, mediators and colliders) and not a validated biomarker panel; the next step is `bio-machine-learning-biomarker-discovery` with independent data.

Done when the top modules are classified as technical or biological candidates and the batch check outcome is stated.
