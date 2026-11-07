## Ligand-Based Pharmacophore (RDKit Pharm3D)

**Goal:** Derive a common pharmacophore from aligned bioactive conformers, then apply that established model to candidate molecules.

**Approach:** Consensus derivation is a separate modeling step: select or generate plausible bioactive conformers, align them using a documented method, identify conserved feature correspondences, and estimate distance bounds or tolerances. RDKit does not provide a single `EmbedPharmacophore` call that performs those steps. `EmbedPharmacophore` instead generates conformations of a molecule that satisfy an already defined pharmacophore.

```python
from rdkit import Chem, Geometry
from rdkit.Chem import ChemicalFeatures
from rdkit.Chem.Pharm3D import EmbedLib, Pharmacophore
from rdkit.RDPaths import RDDataDir
import os

fdef_file = os.path.join(RDDataDir, 'BaseFeatures.fdef')
factory = ChemicalFeatures.BuildFeatureFactory(fdef_file)

# This is an already defined model. Coordinates and bounds must come from a
# validated consensus-derivation workflow or another justified source. RDKit
# requires FreeChemicalFeature objects, not feature-family strings.
# PLACEHOLDER COORDINATES -- do not reuse without deriving from a validated
# workflow. On a real 5-active/5-inactive HIV protease inhibitor retrospective
# set, this literal 2-feature model gave enrichment 0.5 (worse than random).
query_features = [
    ChemicalFeatures.FreeChemicalFeature(
        'Aromatic', Geometry.Point3D(0.0, 0.0, 0.0)),
    ChemicalFeatures.FreeChemicalFeature(
        'Donor', Geometry.Point3D(4.0, 0.0, 0.0)),
]
pharmacophore = Pharmacophore.Pharmacophore(query_features)
pharmacophore.setLowerBound(0, 1, 3.5)
pharmacophore.setUpperBound(0, 1, 5.0)

target = Chem.AddHs(Chem.MolFromSmiles('c1ccc(cc1)CCN'))
can_match, feature_matches = EmbedLib.MatchPharmacophoreToMol(
    target, factory, pharmacophore)
if can_match:
    atom_match = tuple(tuple(matches[0].GetAtomIds())
                       for matches in feature_matches)
    _, embeddings, n_failed = EmbedLib.EmbedPharmacophore(
        target, atom_match, pharmacophore, randomSeed=23, silent=True)
```

`BaseFeatures.fdef` (RDKit-shipped) defines feature SMARTS and is a useful starting feature taxonomy. The code above demonstrates applying an existing two-feature model; it does not infer a consensus model from active compounds.
