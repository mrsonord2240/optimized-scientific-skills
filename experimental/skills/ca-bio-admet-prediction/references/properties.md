## PAINS, BRENK, REOS Filters

ADMET prediction is separate from structural alerts; combine. See `chemoinformatics/substructure-search` for PAINS/BRENK/REOS pattern catalogs.

```python
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams

def alerts(mol, catalogs=('PAINS_A', 'BRENK', 'ZINC')):
    params = FilterCatalogParams()
    for cat in catalogs:
        params.AddCatalog(getattr(FilterCatalogParams.FilterCatalogs, cat))
    catalog = FilterCatalog(params)
    hits = catalog.GetMatches(mol)
    return [h.GetDescription() for h in hits]
```

## Lipinski / Veber / Drug-Likeness

See `chemoinformatics/molecular-descriptors` for full physchem table. Quick filter:

```python
from rdkit.Chem import Descriptors, Lipinski, QED

def druglike_score(mol):
    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    hbd = Lipinski.NumHDonors(mol)
    hba = Lipinski.NumHAcceptors(mol)
    tpsa = Descriptors.TPSA(mol)
    rotbonds = Lipinski.NumRotatableBonds(mol)
    qed = QED.qed(mol)

    lipinski_violations = sum([mw > 500, logp > 5, hbd > 5, hba > 10])
    veber_pass = rotbonds <= 10 and tpsa <= 140
    bbb_simple_screen = tpsa <= 90 and mw <= 500 and hbd <= 3

    return {'MW': mw, 'LogP': logp, 'HBD': hbd, 'HBA': hba,
            'TPSA': tpsa, 'RotBonds': rotbonds, 'QED': qed,
            'Lipinski_violations': lipinski_violations,
            'Veber_pass': veber_pass, 'BBB_simple_screen': bbb_simple_screen}
```
