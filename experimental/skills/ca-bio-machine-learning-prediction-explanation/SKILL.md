---
name: ca-bio-machine-learning-prediction-explanation
category: Data Analysis
description: Use when explaining what an omics classifier relies on with SHAP, LIME or permutation importance, or deciding whether an attribution ranking can be read as biology.
tool_type: python
primary_tool: shap
license: MIT
author: GPTomics
---

# Explaining Omics Classifier Predictions

Pick the one row that matches your question. Read that file and only that file. Run its command before writing any code of your own.
Paths are relative to this Skill's directory.

| Your question | Read |
|---------------|------|
| Which genes drive a tree model (global attribution, correlated genes) | `routes/shap-tree.md` |
| Why one sample was predicted as it was | `routes/lime-local.md` |
| Permutation importance, or the model is linear | `routes/permutation-linear.md` |
| Is the model cheating on batch or a technical signal | `routes/shortcut-audit.md` |
| Which method fits my question (SHAP modes, KernelSHAP, DeepSHAP, LIME) | `routes/choose-method.md` |
| Pick a biomarker panel from the ranking | `bio-machine-learning-biomarker-discovery` |
| Train the model, or check that it generalizes | `ca-bio-machine-learning-omics-classifiers`, `ca-bio-machine-learning-model-validation` |

Rules on every route:

- Input is a CSV with one row per sample, a binary label column and numeric feature columns (loaded by `scripts/omics_data.py`; `--positive` names the positive class, `--drop` removes id and batch columns). Explain only a model that generalizes: the scripts refuse a test AUC under 0.6.
- An attribution describes the model, not biology and not causation. A high-SHAP gene can be a correlate of batch.
- A SHAP or LIME ranking is not a validated biomarker selection: no FDR, no replication, winner's curse.
- Genes co-express in modules and the algorithm decides how to split credit inside one. Aggregate attributions over modules before ranking; the order within a module is not a finding.
- Set the SHAP conditioning explicitly (`feature_perturbation`); the `'auto'` default flips the estimand on whether a background is given.
- State the attribution scale (probability or log-odds) and the background used.

Tested with shap 0.51 (checked environment 0.48-0.51; the script was not rerun), lime 0.2, scikit-learn 1.9.1.
