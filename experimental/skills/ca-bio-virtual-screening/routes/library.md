# Ligand-library screen

```python
import json, re, pandas as pd
from scripts.virtual_screen import virtual_screen
df = pd.read_csv('library.csv')
if not {'name', 'smiles'} <= set(df.columns) or df['name'].duplicated().any():
    raise ValueError('library needs unique name and smiles columns')
if any(not re.fullmatch(r'[A-Za-z0-9_-]+', str(name)) for name in df['name']):
    raise ValueError('names must be safe filename components')
box = json.load(open('box.json'))
virtual_screen('receptor.pdbqt', dict(zip(df['name'], df['smiles'])), box['center'], box['size'], 'screen_output', exhaustiveness=8, seed=42)
```

box.json supplies center and size as three-number arrays. This source helper requires the Vina Python API and writes ligand files, poses and results.csv. It records per-ligand errors; if all affinities are nonnegative it falls back to reporting them, so inspect and flag those modes. Validate the search on controls; see `references/vina.md` and `references/failure-modes.md`.

Done when ranked results, errors, box and settings are retained.
