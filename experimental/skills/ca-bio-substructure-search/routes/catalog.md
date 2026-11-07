# Structural-alert catalogs

```python
import json
from rdkit import Chem
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
params = FilterCatalogParams(); params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_A)
catalog = FilterCatalog(params)
smiles = [x.strip() for x in open('library.smi') if x.strip()]
rows = []
for i, s in enumerate(smiles):
    mol = Chem.MolFromSmiles(s)
    hits = [h.GetDescription() for h in catalog.GetMatches(mol)] if mol else None
    rows.append({'index': i, 'smiles': s, 'alerts': hits})
json.dump(rows, open('alerts.json', 'w'), indent=2)
```

Choose the catalog deliberately; None indicates invalid input. PAINS/BRENK flags require assay-context review, not automatic deletion. See `references/alerts.md`.

Done when alert categories and invalid inputs are retained.
