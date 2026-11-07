# Physicochemical and structural alerts

```python
from rdkit import Chem
from scripts.predict_admet import batch_druglikeness, flag_structural_alerts
molecules = [Chem.MolFromSmiles(x.strip()) for x in open('input.smi') if x.strip()]
batch_druglikeness(molecules).to_csv('properties.csv', index=False)
clean, alerts = flag_structural_alerts(molecules)
import json
json.dump([{'smiles': Chem.MolToSmiles(m), 'alerts': a} for m, a in alerts], open('alerts.json', 'w'), indent=2)
```

Record parse failures separately. PAINS and other alerts flag possible assay interference for review, not automatic exclusions. For thresholds and the simple BBB screen see `references/properties.md`.

Done when properties and alert descriptions are retained for review.
