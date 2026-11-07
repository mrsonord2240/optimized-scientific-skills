## Physicochemical Descriptors

| Descriptor | Source | Range | Drug-like cutoff |
|------------|--------|-------|-------------------|
| MolWt | RDKit `Descriptors.MolWt` | ~50-2000 Da | <=500 (Lipinski) |
| MolLogP (Crippen) | RDKit `Descriptors.MolLogP` | -5 to 8 | <=5 (Lipinski) |
| HBD | `Lipinski.NumHDonors` | 0-10 | <=5 (Lipinski) |
| HBA | `Lipinski.NumHAcceptors` | 0-15 | <=10 (Lipinski) |
| TPSA | `Descriptors.TPSA` (Ertl) | 0-200 A^2 | <=140 (Veber oral); <=90 (BBB+) |
| RotBonds | `Lipinski.NumRotatableBonds` | 0-15 | <=10 (Veber) |
| AromaticRings | `Lipinski.NumAromaticRings` | 0-6 | <=3-4 (Ritchie-Macdonald aromatic ring count) |
| HeavyAtoms | `Descriptors.HeavyAtomCount` | <=50 (lead-like) | |
| FractionCSP3 | `Descriptors.FractionCSP3` | 0-1 | Descriptive; higher sp3 character was associated with clinical progression by Lovering et al. (2009), without a universal cutoff |
| QED | `QED.qed` | 0-1 | Higher is more similar to the reference property distributions; a project may use >=0.5 as a triage heuristic |
| SAscore | `sascorer.calculateScore` (external) | 1-10 | Lower is easier by the model; project cutoffs such as <=4 or >6 require dataset calibration |

**Goal:** Compute a standard physicochemical descriptor panel for drug-likeness filtering and QSAR features.

**Approach:** Combine RDKit `Descriptors`, `Lipinski`, and `QED` calls into a single dict so the caller gets MW, LogP, HBD/HBA, TPSA, rotatable bonds, aromatic rings, fraction sp3, and QED in one pass.

```python
from rdkit.Chem import Descriptors, Lipinski, QED

def physchem(mol):
    return {
        'MolWt': Descriptors.MolWt(mol),
        'MolLogP': Descriptors.MolLogP(mol),
        'HBD': Lipinski.NumHDonors(mol),
        'HBA': Lipinski.NumHAcceptors(mol),
        'TPSA': Descriptors.TPSA(mol),
        'RotBonds': Lipinski.NumRotatableBonds(mol),
        'AromRings': Lipinski.NumAromaticRings(mol),
        'FractionCSP3': Descriptors.FractionCSP3(mol),
        'QED': QED.qed(mol),
    }
```

## Drug-Likeness Rule Sets

| Rule | Constraints | Source |
|------|-------------|--------|
| Lipinski Ro5 | MW<=500, LogP<=5, HBD<=5, HBA<=10 | Lipinski 1997 |
| Veber | RotBonds<=10, TPSA<=140 | Veber 2002 (oral) |
| Ghose | 160<=MW<=480, -0.4<=LogP<=5.6, 40<=MR<=130, 20<=atoms<=70 | Ghose 1999 |
| Egan | LogP<=5.88, TPSA<=131.6 | Egan 2000 |
| Muegge | 200<=MW<=600, -2<=LogP<=5, TPSA<=150, rings<=7, C>4, heteroatoms>1, RotBonds<=15, HBD<=5, HBA<=10 | Muegge 2001 |
| Lead-like | MW<=350, LogP<=3 | Teague 1999 |
| Fragment Ro3 | MW<=300, LogP<=3, HBD<=3, HBA<=3, RotBonds<=3, TPSA<=60 A^2 | Congreve 2003 |
| Pfizer CNS MPO | Six desirability functions: ClogP, ClogD, MW, TPSA, HBD, and pKa | Wager 2010 |

**Use case:** Treat Ro5 and Veber criteria as risk indicators rather than universal hard cutoffs. Doak et al. (2014) analyze orally bioavailable drugs and candidates beyond the Rule of 5, but do not support the claim that approximately 30% of marketed oral drugs violate at least one rule. For CNS prioritization, implement the six-property Wager MPO desirability score rather than replacing it with three hard thresholds.

## QED (Weighted Drug-Likeness)

QED (Bickerton 2012) is a single-number drug-likeness measure (0-1) combining 8 properties (MW, LogP, HBD, HBA, PSA, RotBonds, AromaticRings, structural alerts) via desirability functions.

**Caveat:** QED summarizes desirability functions derived from property distributions of marketed oral drugs; it is not a supervised predictor trained specifically on FDA-approved drugs. It can under-rank large natural-product-like or peptide-like chemistry, while small fragments can be over-ranked because they easily satisfy several desirability functions. Do not use it as the sole filter in either regime; inspect the unfiltered distribution and calibrate against the project objective.
