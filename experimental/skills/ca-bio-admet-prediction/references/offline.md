## Executable Offline Route: ADMET-AI 2.x

`scripts/predict_admet.py` includes an offline batch route. ADMET-AI 2.0.1 ships its model artefacts with the installed package and returns a pandas DataFrame for a list of valid SMILES. It covers common endpoints including hERG, AMES, DILI, Caco-2, BBB, and CYP models; inspect the returned columns in the installed version rather than assuming a fixed endpoint list.

Run the chemistry/input gate before predicting. It rejects malformed, disconnected, inorganic, and metal-containing inputs, and reports maximum Morgan-fingerprint Tanimoto similarity to a **relevant local reference set**. The reference set is a transparent similarity screen, not proof of the model's training applicability domain. Derive its threshold from your reference chemistry with `fit_similarity_threshold`; do not use an arbitrary universal cutoff.

```python
from pathlib import Path
from scripts.predict_admet import (
    assess_applicability_domain,
    fit_similarity_threshold,
    predict_admet_ai,
)

reference_smiles = [line.strip() for line in Path('relevant_reference.smi').read_text().splitlines()
                    if line.strip()]
threshold = fit_similarity_threshold(reference_smiles)  # low-tail, leave-one-out similarity
gate = assess_applicability_domain(smiles_batch, reference_smiles, threshold=threshold)
accepted = gate.loc[gate['decision'] == 'predict', 'smiles'].tolist()

predictions = predict_admet_ai(accepted, endpoints=['hERG', 'AMES', 'DILI', 'BBB_Martins'])
# Point predictions are prioritisation evidence, not safety or regulatory decisions.
```

Inputs marked `manual_review` have organic chemistry but lie outside the local similarity threshold or exceed the configured size limit. Do not silently drop them: retain the gate table and choose experimental data, a more relevant model, or a justified expanded reference set. Inputs marked `reject` are outside this small-molecule organic route and must not be passed to the predictor.

ADMET-AI does not supply a per-prediction uncertainty or formal applicability-domain column in this route. Do not describe its DrugBank percentile columns as uncertainty. If modeled-endpoint uncertainty is required, use a service or an in-house model that actually supplies and documents it.
