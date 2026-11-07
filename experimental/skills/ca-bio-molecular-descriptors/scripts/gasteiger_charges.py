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
