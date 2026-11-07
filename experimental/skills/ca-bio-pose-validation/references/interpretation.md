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

Run `python scripts/pose_qc_batch.py receptor.pdb poses1.sdf poses2.sdf ...` (or import `pose_qc_pipeline(docked_sdfs, receptor_pdb)`); it runs `dock` on each file and keeps the first PB-valid pose per file. For one SDF with strain columns added, use `pose_qc_pipeline` in `scripts/validate_poses.py`.

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
