# Preserve ligand chemistry in pose export

```python
from meeko import PDBQTMolecule, RDKitMolCreate
from rdkit import Chem
pm = PDBQTMolecule.from_file('poses.pdbqt', skip_typing=True)
mol = RDKitMolCreate.from_pdbqt_mol(pm)[0]
with Chem.SDWriter('poses.sdf') as writer:
    for cid in range(mol.GetNumConformers()): writer.write(mol, confId=cid)
```

This reconstruction needs Meeko REMARK SMILES records in the input PDBQT; otherwise retain original molecular identity separately. Open Babel PDBQT conversion can lose bond orders and charge. See `references/errors.md` and `bio-pose-validation`.

Done when exported pose chemistry and physical validity are reviewed.
