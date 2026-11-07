# Partial charges

```python
from rdkit import Chem
from rdkit.Chem import AllChem
```

Use the guarded Gasteiger function in `references/charges.md`; it adds hydrogens, rejects unsupported elements and checks finite charge balance. Match the charge method to the downstream force field.

Done when charges and method limitations are recorded.
