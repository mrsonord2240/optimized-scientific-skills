---
name: bio-shape-similarity
description: Performs 3D shape-based similarity searching using ROCS (OpenEye), USRCAT (ultra-fast), Open3DAlign (RDKit), ESPSim (electrostatic), and ShaEP with explicit handling of Tanimoto-Combo (shape + color), shape vs ECFP4 complementarity, conformer-ensemble searching, alignment optimization, and scaffold hopping. Use when searching for shape-mimicking compounds with different scaffolds, identifying bioisosteric replacements, prospective scaffold hopping, or expanding hit series beyond 2D similarity.
tool_type: python
primary_tool: RDKit
license: MIT
author: GPTomics
---

## Version Compatibility

Reference examples tested with: RDKit 2024.09+ (Open3DAlign and USRCAT); official ShaEP syntax checked against ShaEP 1.4.2; ESPSim example checked against espsim 0.0.1; ROCS/FastROCS/ROCS X are commercial OpenEye products.

Install: `conda install -c conda-forge rdkit`. RDKit provides USRCAT through `rdMolDescriptors`, so no separate `usrcat` package is needed; ShaEP is a separate binary (verify the official current release), ESPSim is `pip install espsim`.

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Shape Similarity

Search for compounds with similar 3D shape (and optionally chemical features) to a query molecule. Shape-based screening complements 2D fingerprint search: it can find scaffold-hopped compounds that ECFP4 misses (different scaffolds with similar shape). ROCS (OpenEye) is the industry-standard commercial tool; Open3DAlign (RDKit), USRCAT (Schreyer & Blundell 2012), and ShaEP are open-source alternatives. Modern best practice combines shape with color (chemical-feature similarity) via Tanimoto-Combo: matches share both shape and pharmacophore feature distribution.

For 2D fingerprint similarity, see `chemoinformatics/similarity-searching`. For pharmacophore search (discrete feature constraints), see `chemoinformatics/pharmacophore-modeling`. For 3D conformer generation, see `chemoinformatics/conformer-generation`.

## Shape Method Taxonomy

| Tool | Speed | Approach | Open-source | Fails when |
|------|-------|----------|-------------|------------|
| ROCS / FastROCS (OpenEye) | Hardware/database/conformer-dependent; vendor reports millions of conformers/s for FastROCS | Gaussian shape + color | No | License and prepared database |
| ROCS X | Trillion-scale reaction/synthon space on Orion | FastROCS plus Bayesian-bandit sampling | No | Commercial cloud workflow |
| USRCAT | Very fast alignment-free descriptor comparison | Moment-based + atom types | Yes | Coarse approximation |
| Open3DAlign (RDKit) | Medium | MMFF atom-type/charge-weighted alignment | Yes | Requires compatible typed 3D structures |
| ShaEP | Benchmark on actual conformers/hardware | Field-based (shape + ESP) | Free binary; inspect license | Requires valid 3D structures and charges for ESP |
| ESPSim | Benchmark on actual workload | Electrostatic + shape | Yes | Limited public benchmarks |
| Phase-Shape (Schrödinger) | commercial | Shape + pharmacophore | No | Commercial |
| USR (original) | Very fast alignment-free comparison | Moment-based only | Yes | No atom-type information |

**Decision:** Select a shape method by matched retrieval/enrichment performance, conformer preparation, throughput, licensing, and score semantics. USRCAT is useful as a fast prefilter; Open3DAlign provides an open alignment method; ROCS/FastROCS provide commercial shape/color workflows.

## Decision Tree by Scenario

| Scenario | Method | Notes |
|----------|--------|-------|
| Large prepared library | USRCAT pre-filter + Open3DAlign rescore | Choose rescore budget from measured retrieval saturation; see `references/usrcat.md`, `references/open3dalign.md` |
| Production VS for scaffold hop | ROCS + color (commercial) | Industry standard |
| Scaffold hopping prospective | Open3DAlign with conformer ensemble | Shape + flexibility; see `references/open3dalign.md` |
| Bioisostere replacement | ROCS color with neutral scoring | Pharmacophore-equivalent matches |
| Patent space carve-out | Shape constraint + 2D dissimilarity | Combine shape + dissimilar scaffold |
| Library diversity assessment | USRCAT k-nearest neighbor | Fast; see `references/usrcat.md` |
| Crystal-bound conformer template | Open3DAlign starting from co-crystal pose | Bioactive shape; see `references/open3dalign.md` |
| Cross-target screening | Shape + pharmacophore feature | Combined screen |
| ESP-relevant pocket, electrostatic bioisosteres | ShaEP or ESPSim | See `references/esp-similarity.md` |

## Tanimoto-Combo Scoring (ROCS Standard)

TanimotoCombo = Tanimoto_shape + Tanimoto_color

- Tanimoto_shape: volume overlap normalized
- Tanimoto_color: pharmacophore feature overlap

Each component is normalized from 0 to 1, so TanimotoCombo ranges from 0 to 2. It is a sum, not an average. Select follow-up thresholds from a relevant benchmark or enrichment study; a single cutoff is not portable across query preparation, color-force-field settings, and library composition.

## Reference Files

Read the file for the method you are running; the decision tree above says which.

| File | Read when |
|------|-----------|
| `references/usrcat.md` | Alignment-free USRCAT descriptors and scoring, the fast pre-filter |
| `references/open3dalign.md` | Open3DAlign (O3A) alignment and scoring, and the conformer-ensemble search (`shape_search_ensemble`) that the scaffold-hop function below calls |
| `references/esp-similarity.md` | Electrostatic-aware shape: ShaEP (mol2 preparation, obabel NaN-coordinate check) and ESPSim |

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

## References

- Hawkins et al., *J. Med. Chem.* 50:74-82 (2007), DOI 10.1021/jm0603365 -- ROCS virtual-screening comparison.
- Schreyer AM, Blundell T. *J. Cheminformatics* 4:27 (2012) -- USRCAT (DOI 10.1186/1758-2946-4-27).
- Vainio, Puranen & Johnson, *J. Chem. Inf. Model.* 49:492-502 (2009), DOI 10.1021/ci800315d -- ShaEP.
- Tosco, Balle & Shiri, *J. Comput. Aided Mol. Des.* 25:777-783 (2011), DOI 10.1007/s10822-011-9462-9 -- Open3DALIGN.
- RDKit O3A API: https://www.rdkit.org/docs/source/rdkit.Chem.rdMolAlign.html
- RDKit shape API: https://www.rdkit.org/docs/source/rdkit.Chem.rdShapeHelpers.html
- ShaEP official documentation/examples: https://cheminformatics.fi/
- OpenEye ROCS X product documentation: https://www.eyesopen.com/rocsx

## Related Skills

- chemoinformatics/molecular-io - Parse query and library
- chemoinformatics/conformer-generation - Generate 3D conformer ensembles
- chemoinformatics/similarity-searching - 2D similarity comparison
- chemoinformatics/pharmacophore-modeling - Pharmacophore alternative
- chemoinformatics/scaffold-analysis - 2D scaffold analysis
- chemoinformatics/virtual-screening - Shape as pre-filter to docking
