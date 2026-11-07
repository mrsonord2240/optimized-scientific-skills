# Standardized QSAR dataset

```python
import json, pandas as pd
from scripts.standardize_library import prepare_qsar_data
rows, statuses = prepare_qsar_data(pd.read_csv('input.csv'), smiles_col='smiles', activity_col='pIC50')
rows.to_csv('standardized.csv', index=False)
json.dump(statuses, open('statuses.json', 'w'), indent=2)
```

Reconcile assay definitions and units before averaging. Default flags replicate activity ranges above 1.0 rather than dropping them. Preserve status counts and original structures; see `references/qsar.md`. For mixed isotopes verify the chosen pipeline's behavior and per-record policy.

Done when deduplicated rows, disagreement flags and status counts are saved.
