# Receptor interaction features

```python
from plip.structure.preparation import PDBComplex
complex = PDBComplex()
complex.load_pdb('complex.pdb')
complex.analyze()
```

Map typed PLIP interactions to ligand features and retain source atom identities. For an apo pocket use the source-documented apo2ph4 workflow; its PML output is not automatically a Pharmit query. See `references/receptor.md` and `references/troubleshooting.md`.

Done when pocket features and their structural evidence are documented.
