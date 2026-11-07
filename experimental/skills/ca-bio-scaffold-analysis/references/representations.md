## Scaffold Representation Taxonomy

| Representation | Origin | Definition | Use case | Fails when |
|----------------|--------|------------|----------|------------|
| Bemis-Murcko scaffold | Bemis & Murcko 1996 | Ring systems + linkers, R-groups stripped | Default chemotype identifier | Linear molecules (no rings) -> empty scaffold |
| Generic framework | Bemis & Murcko 1996 | Bemis-Murcko with all atoms set to C, all bonds single | Topology comparison | Loses heteroatom info; may merge distinct chemotypes |
| Cyclic skeleton (CSK) | Custom RDKit transformation | Ring atoms only, all C, all single | Pure ring-topology view | Loses linker info; not a built-in Murcko option |
| Murcko atom indices | Derived by matching the scaffold to the parent | Parent-molecule atom indices | Programmatic operations | Symmetry can yield multiple equivalent matches |

```python
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold

def all_scaffold_views(smi):
    mol = Chem.MolFromSmiles(smi)
    bm = MurckoScaffold.GetScaffoldForMol(mol)
    bm_smi = Chem.MolToSmiles(bm)

    generic = MurckoScaffold.MakeScaffoldGeneric(bm)
    generic_smi = Chem.MolToSmiles(generic)

    return {
        'bemis_murcko': bm_smi,
        'generic_framework': generic_smi,
    }
```

Example (RDKit 2026.03.6 canonical SMILES): `Cc1ccc(C(=O)NCC2CCCC2)cc1` -> Bemis-Murcko `O=C(NCC1CCCC1)c1ccccc1`; generic `CC(CCC1CCCC1)C1CCCCC1`.

The amount of collapse is library-dependent. For scale only, 1,500 ChEMBL hERG compounds in the reference run produced 857 Bemis-Murcko scaffolds but 583 generic frameworks (1.47-fold fewer groups); its largest generic framework merged four Bemis-Murcko scaffolds. Measure this ratio on the library at hand before using the generic view for a decision.
