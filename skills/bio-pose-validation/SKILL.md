---
name: bio-pose-validation
description: Validates docked / generated protein-ligand poses using PoseBusters physical-validity tests, strain energy quantification, geometric checks (planarity, vdW overlap, bond/angle distortion), and pose-energy reasonableness. Use when QC-ing docking results, comparing classical vs ML docking outputs, or filtering pose lists before SAR analysis.
tool_type: python
primary_tool: PoseBusters
license: MIT
author: GPTomics
---

## Version Compatibility

Install: `pip install posebusters rdkit pandas`

Reference examples tested with: PoseBusters 0.6+, RDKit 2024.09+, pandas 2.2+, posecheck 0.5+ (optional).

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Pose Validation

Test docked or AI-generated protein-ligand poses for physical plausibility. PoseBusters (Buttenschoen et al. 2024) provides geometric, chemical, and energetic checks that flag implausible poses, including non-planar aromatic rings, van der Waals clashes, broken bonds, altered stereochemistry, and unfavorable internal energies. On the Astex Diverse Set, DiffDock achieved 72% RMSD success but only 47% combined RMSD-and-PB-valid success; the size of this gap is dataset- and method-dependent. PB-valid status complements RMSD for downstream SAR, FEP setup, or generative-model training.

For docking, see `chemoinformatics/virtual-screening`. For ML docking specifically, see `chemoinformatics/ml-docking-rescoring`.

## PoseBusters Test Suite

PoseBusters runs ~20 individual checks grouped into:

The thresholds below are the benchmark criteria reported by Buttenschoen et al. (2024). Installed PoseBusters defaults may differ by version and configuration, so record the package version and resolved configuration.

| Check group | What it tests | 2024 benchmark criterion |
|-------------|---------------|-----------|
| Sanity | Ligand chemical sanity | RDKit sanitization passes |
| Bond lengths | Bond lengths within reference | 0.75–1.25 times RDKit distance-geometry bounds |
| Bond angles | 1–3 distances within reference | 0.75–1.25 times RDKit distance-geometry bounds |
| Internal steric | No intra-ligand clash | Pair distance > 0.70 times the RDKit lower bound |
| Aromatic ring planarity | Aromatic rings planar | Maximum deviation from fitted plane <= 0.25 Å |
| Double-bond stereo | Z/E preserved | Match input SMILES |
| Internal energy | Energy relative to generated conformers | UFF energy ratio <= 100 versus the mean of 50 generated, relaxed conformers |
| Volume overlap | vdW overlap with protein | < 7.5% of ligand vdW volume |
| Minimum distance | No severe protein-ligand clash | Distance >= 0.75 times the sum of vdW radii |
| Chirality | R/S preserved from input | Match input SMILES |

Double-bond stereo, chirality, molecular formula, and molecular bond identity are reference checks: they only run when `mol_true` is supplied (`config='redock'`). Without a reference (`config='dock'`), those checks -- and RMSD -- are absent from the results, not merely unreported (verified against PoseBusters 0.6.5: `dock` returns 22 boolean columns, `redock` returns 28, the difference being `rmsd_<=_2å`, `mol_true_loaded`, `molecular_formula`, `molecular_bonds`, `double_bond_stereochemistry`, and `tetrahedral_chirality`).

A pose passing ALL tests is "PB-valid". Combined PB-valid + RMSD <= 2 Å is the modern criterion.

## When to Apply PoseBusters

| Workflow | PoseBusters use | Action |
|----------|-----------------|--------|
| Self-docking (validating method) | Required | Compare PB-valid + RMSD <= 2A |
| Cross-docking | Required | PB-valid + RMSD <= 2A; account for protein flexibility |
| Virtual screening top hits | Required | Filter to PB-valid before MM/GBSA / FEP |
| AI docking (DiffDock, etc.) | Required for a fair benchmark | Report the dataset-specific PB-valid and combined success rates |
| Generated ligand poses | Recommended | Measure chemical and geometric validity rather than assuming it |
| Boltz-2 / AlphaFold3 ligand poses | Recommended | Benchmark validity on the relevant complexes; do not infer a failure frequency from DiffDock |
| Production FEP setup | Required | Inspect pose validity and ligand strain before system preparation |

## PoseBusters Usage

```python
from posebusters import PoseBusters

bust = PoseBusters(config='redock')

results = bust.bust(
    mol_pred='predicted.sdf',
    mol_true='reference.sdf',
    mol_cond='receptor.pdb',
)
```

Common configurations and their included checks are:

| Config | Includes | When to use |
|--------|----------|-------------|
| `redock` | All checks + RMSD vs reference + protein vdW overlap | Self-docking benchmarks, retrospective validation |
| `dock` | All checks that need no reference; drops RMSD and the reference-dependent checks listed under the check table | Blind docking, prospective virtual screening |
| `mol` | Intra-ligand only (sanitization, bond lengths/angles, internal clash, ring and double-bond flatness, energy); no stereo or chirality checks, which need `mol_true` | Conformer QC; no protein context |

PoseBusters also ships additional and faster configurations in some releases. Treat the table as a workflow guide, not an exhaustive registry, and inspect the configurations available in the installed version.

Output: a DataFrame with one row per pose, metadata columns, and boolean pass/fail columns for the checks enabled by the selected configuration. Reference-dependent fields such as RMSD and the exact check-column names vary by configuration and version; inspect `results.columns` rather than relying on a fixed exhaustive list.

## Python Library API

**Goal:** Programmatically validate a docked-pose SDF against a receptor PDB and produce a PB-valid filter.

**Approach:** Instantiate `PoseBusters(config='dock')`, call `bust()` on the SDF + PDB pair, and AND-aggregate all boolean check columns into a single `pb_valid` flag.

Run `run_posebusters(pred_sdf, receptor_pdb)` from `examples/validate_poses.py` (it uses `dock`, or `redock` when a reference SDF is given, and adds the `pb_valid` column). Check columns are the boolean ones excluding `rmsd*`.

## Strain Energy Quantification

Beyond binary PB-valid, quantitative strain energy distinguishes "marginal" from "egregious" poses.

**Goal:** Quantify how far each docked pose is from its lowest-energy free conformer in MMFF94 energy units.

**Approach:** Generate a reference conformer ensemble (ETKDGv3 + MMFF94), make the docked and reference molecules chemically consistent by adding explicit hydrogens to both, relax only the added docked-pose hydrogens while fixing all heavy atoms, take the lowest sampled reference energy as baseline, and report `docked_energy - min_ref_energy` as a relative strain diagnostic. This is not a rigorous solution-phase conformational free energy.

Run `ligand_strain_mmff(docked_sdf, n_ref_conf=20)` from `examples/validate_poses.py`; it returns a DataFrame with `pose_idx`, `strain_kcal` and a `note` for poses it could not score.

Interpret relative MMFF strain in the context of ligand chemistry, conformer-sampling coverage, and force-field support. Boström et al. (1998) found a conformational energy penalty of no more than 3 kcal/mol for about 70% of 33 protein-bound ligands; that result does not establish a universal acceptance cutoff. Treat unusually high values as a prompt for inspection or use a project-defined threshold validated for the series.

## vdW Overlap with Protein

The 2024 benchmark criterion limits protein-ligand overlap to 7.5% of the ligand vdW volume, using protein radii scaled by 0.8. PoseBusters' `bust(...)` computes this check; do not substitute an unvalidated pairwise-distance sketch for its volume calculation.

## Aromatic Ring Planarity

Run `python scripts/aromatic_planarity.py poses.sdf` (or import `aromatic_planarity(mol)` from it); it prints the maximum out-of-plane deviation per pose.

This reimplementation approximates PoseBusters' internal flatness computation; it does not reproduce it exactly. Verified against installed PoseBusters 0.6.5 (`aromatic_ring_flatness`) on a displacement series of the same pose: this formula's deviation reliably passes at 0.29 Å and reliably fails at 0.45 Å and above -- roughly double the 0.25 Å figure sometimes quoted for this check. Treat `bust()`'s own `aromatic_ring_flatness` column as authoritative for pass/fail; use this snippet only as a supplementary diagnostic, not a gate.

## Model-Specific Failure Diagnosis

Do not assign a mechanism from the model name or a failed PoseBusters column alone. For DiffDock-L, EquiBind, TANKBind, Boltz, AlphaFold3, or another pose generator, report the observed failed checks on the evaluated dataset, inspect the structures, and compare against the method's documented constraints. A chirality, planarity, bond-geometry, or clash failure may justify filtering or a validated constrained-relaxation protocol, but relaxation must be checked for displacement of the binding mode.

### High strain after Vina docking

**Trigger:** Highly constrained pocket; flexible ligand.

**Symptom:** Relative strain is an outlier for the chemical series even though the pose passes the enabled geometric checks.

**Fix:** Inspect conformer-sampling coverage and force-field support. Compare additional docking or constrained-relaxation settings under a project-validated protocol rather than applying a universal strain or exhaustiveness cutoff.

## Reconciliation: PoseBusters vs RMSD

| RMSD <= 2A | PB-valid | Action |
|------------|----------|--------|
| Yes | Yes | Physically plausible and close to the reference; still validate suitability for the downstream task |
| Yes | No | Close to the reference but fails an enabled plausibility check; inspect the failure and any validated relaxation |
| No | Yes | Physically plausible but different from the reference; investigate alignment, protein state, and alternative binding modes |
| No | No | Different from the reference and fails an enabled plausibility check; inspect both causes before deciding whether to reject |

On the Astex Diverse Set reported by Buttenschoen et al. (2024), DiffDock's top-pose success fell from 72% by RMSD <= 2 Å alone to 47% when PB-validity was also required: a 25-percentage-point gap. Do not generalize that result to a fixed failure rate on other datasets.

## Integration into VS Pipeline

Run `python scripts/pose_qc_batch.py receptor.pdb poses1.sdf poses2.sdf ...` (or import `pose_qc_pipeline(docked_sdfs, receptor_pdb)`); it runs `dock` on each file and keeps the first PB-valid pose per file. For one SDF with strain columns added, use `pose_qc_pipeline` in `examples/validate_poses.py`.

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| Rows or expected checks are missing | Input loading failed or the selected configuration omits those checks | Inspect the returned DataFrame, loading-status columns, input format, and installed configuration |
| RMSD not computed | No reference provided | Pass `mol_true` parameter |
| All checks pass for invalid pose | Wrong receptor file format | Use PDB with hydrogens; PDBQT may not work |
| vdW overlap false positive on covalent | Covalent bond counted as clash | Use covalent docking-specific validation |
| Strain calculation slow | Too many reference conformers | Reduce `n_ref_conf` to 5-10 |
| PoseBusters config error | Wrong or version-incompatible config name | Inspect the installed configuration registry; `redock`, `dock`, and `mol` are common configurations |
| posecheck unavailable | Different tool, similar purpose | `pip install posecheck` for alternative |
| Otherwise-reasonable pose rejected on one borderline check | Binary PB-valid used as an absolute reject | Use PoseBusters as a filter, not an absolute reject; inspect which check failed and by how much before discarding the pose |

## References

- Buttenschoen M, Morris GM, Deane CM. "PoseBusters: AI-based docking methods fail to generate physically valid poses or generalise to novel sequences." *Chem. Sci.* 15:3130–3139 (2024). DOI: 10.1039/D3SC04185A.
- Boström J, Norrby PO, Liljefors T. "Conformational energy penalties of protein-bound ligands." *J. Comput.-Aided Mol. Des.* 12:383–396 (1998). DOI: 10.1023/A:1008007507641.
- Corso G et al. "DiffDock: Diffusion Steps, Twists, and Turns for Molecular Docking." *ICLR* (2023). OpenReview: https://openreview.net/forum?id=kKF8_K-mBbS.
- Stärk H et al. "EquiBind: Geometric Deep Learning for Drug Binding Structure Prediction." *PMLR* 162:20503–20521 (2022). https://proceedings.mlr.press/v162/stark22b.html.
- Lu W et al. "TankBind: Trigonometry-Aware Neural NetworKs for Drug-Protein Binding Structure Prediction." *NeurIPS* 35 (2022). Official repository: https://github.com/luwei0917/TankBind.
- Abramson J et al. "Accurate structure prediction of biomolecular interactions with AlphaFold 3." *Nature* 630:493–500 (2024). DOI: 10.1038/s41586-024-07487-w.
- Boltz official repository and documentation: https://github.com/jwohlwend/boltz.
- PoseBusters documentation, Python API: https://posebusters.readthedocs.io/en/latest/api.html.

## Related Skills

- chemoinformatics/virtual-screening - Source of poses to validate
- chemoinformatics/ml-docking-rescoring - DiffDock, EquiBind, TANKBind validation
- chemoinformatics/molecular-io - SDF format handling
- chemoinformatics/conformer-generation - Generate reference conformer ensemble for strain
- chemoinformatics/free-energy-calculations - PoseBusters-valid poses for FEP input
- chemoinformatics/covalent-design - Covalent pose validation
