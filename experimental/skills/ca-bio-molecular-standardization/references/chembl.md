## ChEMBL Structure Pipeline (Reference Implementation)

ChEMBL's standardization is the most widely-used reference. The Python package `chembl_structure_pipeline` exposes the validated pipeline.

**Goal:** Apply the industry-reference ChEMBL standardization pipeline to a SMILES.

**Approach:** Parse SMILES with RDKit, run `standardize_mol` (sanitize, normalize, and standardize charges), then `get_parent_mol` (strip salts/counter-ions), and emit canonical SMILES. Add `rdMolStandardize.TautomerEnumerator` separately only when the project requires tautomer canonicalization.

```python
from chembl_structure_pipeline import standardize_mol, get_parent_mol
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize

def chembl_pipeline(smi, multi_fragment_policy='flag'):
    try:
        mol = Chem.MolFromSmiles(smi)
    except (TypeError, ValueError):
        return None, 'parse_failure'
    if mol is None:
        return None, 'parse_failure'
    standardized = standardize_mol(mol)
    parent, exclude = get_parent_mol(standardized)
    if exclude:
        return None, 'excluded_by_chembl'
    if len(Chem.GetMolFrags(parent)) > 1:
        if multi_fragment_policy == 'flag':
            return None, 'multi_fragment_parent'
        if multi_fragment_policy == 'largest_fragment':
            parent = rdMolStandardize.LargestFragmentChooser(
                preferOrganic=True
            ).choose(parent)
            return Chem.MolToSmiles(parent), 'ok_largest_fragment_fallback'
        raise ValueError("multi_fragment_policy must be 'flag' or 'largest_fragment'")
    return Chem.MolToSmiles(parent), 'ok'
```

**`standardize_mol`:** sanitize, normalize functional groups, and standardize charges; returns one RDKit molecule.

**`get_parent_mol`:** strip salts/counter-ions and choose the parent; returns `(parent_mol, exclude_flag)`.

Output: canonical SMILES of the selected parent after the ChEMBL transformations, or an explicit status when parsing, exclusion, or parent selection prevents a safe result. `get_parent_mol` can return more than one fragment when every input fragment is on ChEMBL's salt list (for example sodium acetate or choline chloride). The default `multi_fragment_policy='flag'` rejects this ambiguous, unstripped parent; use the explicit RDKit `largest_fragment` fallback only when its loss of co-crystal information is acceptable. Neutralizable acid/base sites may be normalized, but permanent or otherwise non-removable charges can remain; do not assume every emitted parent is neutral.
