## Open3DAlign (RDKit)

Open3DAlign uses MMFF atom types and partial charges to find an atom-based 3D alignment:

**Goal:** Align a target molecule onto a query in 3D and score volume overlap with Open3DAlign.

**Approach:** Build 3D structures for query and target, run `GetO3A`, and call `Align()` to transform the probe in place. `Score()` is the unnormalized O3A objective, not a shape Tanimoto or ROCS TanimotoCombo. If a normalized shape similarity is required, compute `1 - rdShapeHelpers.ShapeTanimotoDist(...)` after alignment.

```python
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign, rdShapeHelpers

query = Chem.MolFromSmiles('CCC(=O)Nc1ccccc1')
query = Chem.AddHs(query)
AllChem.EmbedMolecule(query, AllChem.ETKDGv3())

target = Chem.MolFromSmiles('CCC(=O)Nc1ccc(F)cc1')
target = Chem.AddHs(target)
AllChem.EmbedMolecule(target, AllChem.ETKDGv3())

O3A = rdMolAlign.GetO3A(target, query)
rmsd = O3A.Align()  # aligns target to query in place
o3a_score = O3A.Score()
shape_tanimoto = 1.0 - rdShapeHelpers.ShapeTanimotoDist(target, query)
```

`GetO3A` finds an alignment between conformers; `Align()` applies it and returns RMSD. Keep `o3a_score` and normalized `shape_tanimoto` distinct in outputs. `Align()` mutates the probe's coordinates: copy the probe (`Chem.Mol(target)`) or store the best-aligned conformer if the coordinates are part of the output.

**Open3DAlign vs ROCS:** Open3DAlign is open-source and competitive on small benchmarks; slower than ROCS at scale.

## Conformer-Ensemble Shape Searching

For each library molecule, generate ensemble of conformers; pick best-shape conformer:

**Goal:** Run shape-similarity search over a conformer ensemble per library molecule so bound-conformer-like shapes are recovered.

**Approach:** For each library molecule, reject disconnected fragments (salts have no single shape to compare), add hydrogens, embed n_conf conformers with ETKDGv3, MMFF-optimize, drop only the conformers that fail to converge (not the whole molecule), score the surviving conformers against the query with Open3DAlign, and keep the best score per molecule. Report every molecule that is dropped and why -- do not let a molecule silently disappear from the results.

```bash
python scripts/shape_search_ensemble.py --query 'CC(=O)Nc1ccc(C(=O)c2ccccc2)cc1' --library lib.smi --n-conf 20
```

Prints `SMILES<TAB>shape=<best Tanimoto>` best-first, plus a WARNING for every molecule that lost conformers and a list of every dropped molecule with its reason. The same functions import as `from shape_search_ensemble import shape_search_ensemble` (returns `(mol, best_shape)` pairs; `seed=42` fixes the embedding).

**Critical:** Results depend on conformer coverage. Use an ensemble sized and validated for the library and query (20 conformers is a repository starting budget, not a universal minimum) rather than assuming one conformer is representative. A molecule that loses some but not all of its conformers to non-convergence is still scored on the survivors; only a molecule with zero usable conformers is dropped, and every drop is printed with its SMILES and reason -- never assume "fewer hits than input molecules" means the missing ones failed cleanly.
