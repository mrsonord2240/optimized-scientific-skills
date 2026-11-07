## ADMETlab 3.0 API (optional hosted route)

ADMETlab 3.0 reports 119 platform features: 77 prediction models, 34 computed physicochemical properties, and 8 medicinal-chemistry rules. The modeled endpoints include prediction uncertainty; do not imply that computed properties and rules have model uncertainty.

**Goal:** Obtain the platform's 119 features for a batch of SMILES, including uncertainty for the 77 modeled endpoints, using the hosted service.

**Approach:** Follow the live ADMETlab 3.0 API tutorial to wash molecules, submit batch predictions, and retrieve the returned results. The 2024 paper documents API/batch support and modeled-endpoint uncertainty. Obtain current rate limits, routes, payloads, task identifiers, and output contracts from the live official documentation rather than attributing them to the paper or hard-coding an unofficial example.

```python
from scripts.predict_admet import load_admetlab_results

# After submitting with the current official API example, validate its CSV.
# Supply the current task-ID column name published by that API contract.
results = load_admetlab_results(
    'admetlab3_results.csv', task_id_column='task_id'
)
# Preserve the validated uncertainty columns and task identifier downstream.
```

The package deliberately does not bundle a request URL, payload, rate limit, or task identifier because those must come from the current official service documentation. This is an optional hosted integration, not the Skill's executable prediction route. If the service's current output cannot meet the loader contract, stop and reconcile the contract instead of treating an unrelated CSV as predictions.

ADMETlab 3.0 endpoints (sample):
- Absorption: Caco-2 permeability (logPapp), HIA (%), Pgp inhibitor/substrate, MDCK
- Distribution: BBB+, PPB (%), VDss (L/kg), Fu (fraction unbound)
- Metabolism: CYP1A2/2C9/2C19/2D6/3A4 inhibitor / substrate
- Excretion: CL (mL/min/kg), T1/2 (h)
- Toxicity: hERG, AMES, hepatotoxicity (DILI), carcinogenicity, immunotoxicity, mutagenicity, respiratory, skin, eye, cardiotoxicity, mitochondrial, NR-AR, NR-ER, SR-MMP
- Drug-likeness: Lipinski, Veber, Ghose, Egan, Muegge, QED, SAscore
