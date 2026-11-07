## Open Babel for MOL2 / PDBQT

RDKit's MOL2 parser is incomplete (SYBYL atom-type sets differ). Open Babel is more robust for MOL2 and PDBQT.

```python
from openbabel import pybel
from rdkit import Chem
from rdkit.Chem.inchi import MolToInchi

mols = list(pybel.readfile('mol2', 'ligands.mol2'))
for mol in mols:
    smi = mol.write('smi').strip().split()[0]
    rdkit_mol = Chem.MolFromSmiles(smi)
    if rdkit_mol is None:
        raise ValueError(f'Open Babel produced an unreadable SMILES: {smi!r}')
    inchi = MolToInchi(rdkit_mol)
```

Use Open Babel for MOL2/PDBQT conversion, but generate registry InChI with RDKit as above. Open Babel's InChI support is build-dependent; if an Open Babel InChI route is required, verify it first with `obabel -L formats` and handle an absent `inchi` format rather than assuming `mol.write('inchi')` works.

For docking output PDBQT, use Open Babel rather than RDKit:
```python
import subprocess
subprocess.run(['obabel', 'docked.pdbqt', '-O', 'docked.sdf'], check=True)
```

## InChI for Canonical Identity

InChI is a standardized, canonical structure identifier designed for cross-database and cross-toolkit interoperability (Heller et al. 2015; O'Boyle 2012). Standard InChI normalizes many mobile-hydrogen tautomers and has limitations for some metal stereochemistry; non-standard options such as `/FixedH` can distinguish additional representations. InChIKey is a fixed-length hash, so use full InChI or standardized structures when a suspected collision must be resolved (InChI Trust technical FAQ).

```python
from rdkit.Chem.inchi import MolToInchi, MolToInchiKey, InchiToInchiKey

mol = Chem.MolFromSmiles('c1ccc2c(c1)cccc2')
inchi = MolToInchi(mol)
key = MolToInchiKey(mol)

inchi_fixedH, aux_info = Chem.MolToInchiAndAuxInfo(mol, options='/FixedH')
```

**Caveat:** Two molecules with identical std InChI may be different tautomers. Use `/FixedH` for tautomer-distinguishing InChI when needed.
