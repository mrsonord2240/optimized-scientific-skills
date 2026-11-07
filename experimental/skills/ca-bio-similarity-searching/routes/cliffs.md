# Activity cliffs

```python
import pandas as pd
from scripts.activity_cliffs import activity_cliffs
df = pd.read_csv('activity.csv')
pairs = activity_cliffs(df, sim_threshold=0.85, activity_gap=2.0, activity_col='pIC50', include_chirality=True)
pd.DataFrame(pairs, columns=['row_i','row_j','similarity','activity_gap']).to_csv('cliffs.csv', index=False)
```

The helper is `scripts/activity_cliffs.py`. Thresholds are starting settings; standardize salts and assay units before comparing. Examine measurement uncertainty and stereochemistry; see `references/mcs-cliffs.md` and `references/errors.md`.

Done when pair indices, settings and limitations are saved.
