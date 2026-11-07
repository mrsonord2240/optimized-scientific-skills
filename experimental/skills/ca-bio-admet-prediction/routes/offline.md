# Offline ADMET prediction

```python
from scripts.predict_admet import fit_similarity_threshold, assess_applicability_domain, predict_admet_ai
reference = [x.strip() for x in open('relevant_reference.smi') if x.strip()]
inputs = [x.strip() for x in open('input.smi') if x.strip()]
gate = assess_applicability_domain(inputs, reference, threshold=fit_similarity_threshold(reference))
gate.to_csv('gate.csv', index=False)
accepted = gate.loc[gate.decision == 'predict', 'smiles'].tolist()
if accepted:
    predict_admet_ai(accepted).to_csv('predictions.csv')
```

The helper is `scripts/predict_admet.py`. The reference-similarity gate is not the model's formal training domain. Keep rejected and manual-review inputs. Inspect installed endpoint names and version; the 2.x output is not interchangeable with v1. See `references/offline.md` and `references/selection.md`.

Done when the gate and any predictions are saved with model and reference provenance.
