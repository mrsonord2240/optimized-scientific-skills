# Endpoint interpretation

```python
from scripts.predict_admet import fit_similarity_threshold, assess_applicability_domain, predict_admet_ai
reference = [x.strip() for x in open('relevant_reference.smi') if x.strip()]
inputs = [x.strip() for x in open('input.smi') if x.strip()]
gate = assess_applicability_domain(inputs, reference, threshold=fit_similarity_threshold(reference))
gate.to_csv('gate.csv', index=False)
accepted = gate.loc[gate.decision == 'predict', 'smiles'].tolist()
if accepted:
    predict_admet_ai(accepted).to_csv('endpoint_predictions.csv')
```

The helper is `scripts/predict_admet.py`; inspect available endpoint columns instead of assuming task names. Retain rejected and manual-review records. Review hERG and CYP assay definitions and exposure in `references/endpoints.md`; for disagreements and failure modes use `references/troubleshooting.md`. Point predictions do not supply calibrated uncertainty.

Done when the gate and endpoint results are saved with applicability and uncertainty limitations.
