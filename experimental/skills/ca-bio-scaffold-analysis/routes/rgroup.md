# R-group decomposition

```python
import json, pandas as pd
from scripts.rgroup import decompose_series
df = pd.read_csv('library.csv')
core = open('core.smi').readline().strip()
rows, matched, unmatched = decompose_series(df['smiles'].tolist(), core)
json.dump({'rows': rows, 'matched_positions': matched, 'unmatched_positions': unmatched}, open('rgroups.json', 'w'), indent=2)
```

The helper is `scripts/rgroup.py`. Use labeled attachment points and inspect symmetry-dependent assignments and unmatched rows; see `references/rgroup.md`.

Done when R-group assignments retain their input positions and unmatched indices.
