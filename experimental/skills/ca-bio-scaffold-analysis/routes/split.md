# Scaffold-isolated split

```python
import json, pandas as pd
from scripts.scaffold_split import scaffold_split
df = pd.read_csv('library.csv')
train, test, diagnostics = scaffold_split(df, smiles_col='smiles', train_frac=0.8, seed=42, return_diagnostics=True)
train.to_csv('train.csv', index=False)
test.to_csv('test.csv', index=False)
json.dump(diagnostics, open('split_diagnostics.json', 'w'), indent=2)
```

Whole scaffold groups remain together; check achieved sizes, overlap and endpoint balance. Scaffold split is one deployment question, not a substitute for time or prospective validation. See `references/clustering.md`.

Done when both partitions and diagnostics are saved.
