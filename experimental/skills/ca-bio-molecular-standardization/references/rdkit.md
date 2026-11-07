## Full Standardization with rdMolStandardize

For more granular control or non-ChEMBL workflows.

**Goal:** Execute each standardization step explicitly to control salt stripping, charge handling, tautomer canonicalization, and isotope normalization.

**Approach:** Run the 8-stage pipeline (sanitize, largest fragment, normalize, uncharge, tautomer canonicalize, isotope strip, stereo standardize, canonical SMILES) sequentially with `rdMolStandardize` primitives.

```python
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize

def full_standardize(smi, keep_isotopes=False):
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None

    Chem.SanitizeMol(mol)

    largest = rdMolStandardize.LargestFragmentChooser(preferOrganic=True)
    mol = largest.choose(mol)

    normalizer = rdMolStandardize.Normalizer()
    mol = normalizer.normalize(mol)

    uncharger = rdMolStandardize.Uncharger(canonicalOrder=True)
    mol = uncharger.uncharge(mol)

    enumerator = rdMolStandardize.TautomerEnumerator()
    mol = enumerator.Canonicalize(mol)

    if not keep_isotopes:
        for atom in mol.GetAtoms():
            atom.SetIsotope(0)

    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    return Chem.MolToSmiles(mol)
```

**`canonicalOrder=True`** makes the uncharger choose neutralization sites in canonical order when more than one equivalent site is available. It does not itself decide whether a permanent charge is retained; inspect charge-sensitive structures and keep `force=False` unless a documented policy requires otherwise.

For a mixed tracer library, apply the isotope choice per row rather than forcing one setting across every compound. Store the source decision in a boolean column and pass it to `full_standardize`:

```python
df['standardized_smiles'] = [
    full_standardize(smi, keep_isotopes=keep)
    for smi, keep in zip(df['smiles'], df['keep_isotopes'])
]
```

Require `keep_isotopes` to be boolean, preserve it in the output, and record the policy with the exported dataset. Use a library-wide setting only when the entire collection has the same tracer policy.
