## Shape vs ECFP4 Complementarity

| Shape result | ECFP4 result | Interpretation |
|--------------|--------------|----------------|
| High | High | Close analog candidate |
| High | Low | Scaffold-hop candidate |
| Low | High | Similar 2D chemotype in a different sampled shape |
| Low | Low | Unrelated by these representations |

Calibrate “high” and “low” on a task-relevant reference set; do not treat the illustrative defaults below as universal scientific cutoffs.

The shape >> ECFP4 quadrant is the scaffold-hopping gold:

**Goal:** Identify scaffold-hop candidates that are 3D-shape-similar but 2D-chemotype-dissimilar to the query.

**Approach:** Run the conformer-ensemble shape search (`shape_search_ensemble`, `references/open3dalign.md`), keep hits above a shape Tanimoto cutoff, then retain only those whose ECFP4 Tanimoto to the query is below an ECFP4 dissimilarity cutoff.

```bash
python scripts/shape_search_ensemble.py --query 'CC(=O)Nc1ccc(C(=O)c2ccccc2)cc1' --library lib.smi --scaffold-hop --shape-threshold 0.7 --ecfp-threshold 0.5
```

The 0.7 / 0.5 cutoffs are repository starting defaults only; calibrate both on a task-relevant active/decoy or retrieval benchmark before making decisions. ECFP4 here is Morgan radius 2, 2048 bits (compare only fingerprints built with identical settings). Import `scaffold_hop_candidates` from the script for Python use.

## Per-Tool Failure Modes

### Disconnected-fragment (salt) input -- silent success on meaningless shapes

**Trigger:** Library SMILES with more than one component, e.g. a salt (`[Fe+2].[Cl-].[Cl-]`) or a counter-ion-carrying export.

**Mechanism:** `AllChem.MMFFHasAllMoleculeParams(mol)` returns `True` for disconnected fragments -- unbonded single-atom pieces have nothing for MMFF to fail to parameterize -- so embedding and MMFF optimization both "succeed" and produce conformers with no real combined shape to compare.

**Symptom:** No exception anywhere; the molecule gets a shape score as if it were a normal structure.

**Fix:** Check `len(Chem.GetMolFrags(mol)) > 1` before embedding and reject or warn on multi-fragment input; `shape_search_ensemble` (`references/open3dalign.md`) does this. Salt-strip upstream (see `chemoinformatics/molecular-standardization`) if the intent is to compare the parent structure's shape.

### USRCAT -- false positive on small molecules

**Trigger:** Library has many fragment-sized compounds.

**Mechanism:** USRCAT moments dominated by overall shape; small molecules look "similar" if shape resemble.

**Symptom:** Many fragment hits; not pharmacophore-relevant.

**Fix:** Calibrate size/property filters on the retrieval task and rescore selected hits with an alignment or feature-aware method.

### Open3DAlign -- slow on large library

**Trigger:** Million-compound library, full alignment.

**Mechanism:** Open3DAlign is iterative; O(N) per molecule.

**Symptom:** Hours of compute.

**Fix:** Pre-filter with USRCAT and choose the rescore budget from measured throughput and retrieval saturation.

### Shape only -- wrong stereochemistry match

**Trigger:** Mirror-image of correct binder.

**Mechanism:** Shape-only scoring may insufficiently penalize stereochemical alternatives even though a rigid rotational overlay is not generally invariant to mirror reflection.

**Symptom:** Enantiomer of inactive scores as hit.

**Fix:** Validate hits by 3D pose; check stereochemistry.

### ROCS color -- bioisostere missed

**Trigger:** -COOH replaced by -SO3H or tetrazole.

**Mechanism:** Default color types may not equate these bioisosteres.

**Symptom:** Known bioisostere doesn't score high.

**Fix:** Validate the color-force-field treatment for the bioisostere and compare shape, color, and pharmacophore evidence separately.

### Conformer not bioactive

**Trigger:** Library compound generated conformer is not the bound conformation.

**Mechanism:** ETKDGv3 generates plausible conformers; bound conformer may be higher energy.

**Symptom:** Known active doesn't shape-match query.

**Fix:** Use larger conformer ensemble; weight by Boltzmann; or use CREST + GFN2-xTB for high-quality sampling.

### Field-based methods slower

**Trigger:** ShaEP or ESPSim on production library.

**Mechanism:** Field-based methods compute Gaussian fields per molecule.

**Symptom:** Field calculation or alignment dominates runtime on the prepared library.

**Fix:** Use as second-stage rescore; not primary screen.

## Reconciliation: Shape vs Pharmacophore

| Aspect | Shape | Pharmacophore |
|--------|-------|----------------|
| Representation | Volume distribution | Discrete features in space |
| Captures | Overall bulk | Interaction-relevant features |
| Speed | Fast (USRCAT) to medium (Open3DAlign) | Fast |
| Specificity | Task- and query-dependent | Task- and feature-definition-dependent |
| False positive rate | Measure on a matched benchmark | Measure on a matched benchmark |
| Best for | Scaffold hopping initial | Scaffold hopping refinement |

Shape and pharmacophore searches make different approximations. Compare them alone and in sequence on a matched active/decoy or retrieval benchmark before assigning recall/precision roles.

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| Open3DAlign RMSD is near 0 | Near-exact O3A alignment | Treat as a successful alignment; evaluate the O3A and shape scores separately |
| USRCAT vector all zeros | Mol has no 3D coords | Generate conformer first |
| Shape Tanimoto > 1 | Raw O3A or TanimotoCombo mislabeled as shape Tanimoto | Shape Tanimoto is 0-1; O3A is unnormalized and ROCS TanimotoCombo is 0-2 |
| ROCS very slow | Sequential processing | Use parallel batching |
| Shape match but no docking pose | Wrong binding pose | Use docking on top shape hits, not shape alone |
| Missing co-crystal template | Apo or AlphaFold-only structure | Use ligand-based pharmacophore + shape |
| ShaEP returns no hits | Strict tolerance | Loosen overlap thresholds |
