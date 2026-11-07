## Partial Charge Methods

| Method | Software | Cost | Accuracy | Use for |
|--------|----------|------|----------|---------|
| Gasteiger-Marsili | RDKit, Open Babel | Fast | Empirical, rough | Charge-aware preparation or models that explicitly require Gasteiger charges; Vina/Vinardo scoring itself does not require assigned atom charges |
| MMFF94 | RDKit | 0.1s/mol | Force-field consistent | MMFF energy, conformer ranking |
| AM1-BCC | antechamber (AmberTools) | ~10s/mol | Semi-empirical | MD setup, FEP, GAFF |
| RESP | psi4, Gaussian | minutes/mol | Restrained fit to a quantum-mechanical ESP; protocol-specific | Force-field workflows parameterized for that RESP protocol |
| OpenFF Recharge | openff-recharge | Workflow-dependent | Framework for generating/retrieving QC ESP data and fitting library charges, BCCs, RESP charges, or virtual sites | Developing or evaluating charge models; it is not one charge-assignment method |

```python
import math
from rdkit import Chem
from rdkit.Chem import AllChem

# Conservative organic-subset check: do not interpret a finite value for a
# metal or other unsupported element as a validated Gasteiger parameter.
GASTEIGER_ORGANIC_ELEMENTS = {'H', 'C', 'N', 'O', 'F', 'P', 'S', 'Cl', 'Br', 'I'}

def gasteiger_charges(mol):
    unsupported = sorted({atom.GetSymbol() for atom in mol.GetAtoms()
                          if atom.GetSymbol() not in GASTEIGER_ORGANIC_ELEMENTS})
    if unsupported:
        raise ValueError(
            'Gasteiger charges are not validated here for: ' + ', '.join(unsupported))
    mol_h = Chem.AddHs(Chem.Mol(mol))
    AllChem.ComputeGasteigerCharges(mol_h)
    charges = [float(atom.GetProp('_GasteigerCharge')) for atom in mol_h.GetAtoms()]
    if not all(math.isfinite(charge) for charge in charges):
        raise ValueError('Gasteiger charge calculation returned a non-finite value')
    if not math.isclose(sum(charges), Chem.GetFormalCharge(mol_h), abs_tol=1e-4):
        raise ValueError('Gasteiger charges did not balance to the formal charge')
    return mol_h, charges
```

`_GasteigerCharge` on heavy atoms alone does not balance when hydrogens are implicit: their charge is stored in `_GasteigerHCharge`. The helper makes hydrogens explicit and sums every atom's `_GasteigerCharge` instead. RDKit can return a finite formal-charge-like number for metals (for example, iron) without signalling parameter coverage, so reject non-organic elements rather than treating that number as a usable partial charge.

**Critical:** Charge method must match downstream. Gasteiger charges in an AMBER MD run violate the assumptions of the protein force field.
