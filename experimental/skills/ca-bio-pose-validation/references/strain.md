## Strain Energy Quantification

Beyond binary PB-valid, quantitative strain energy distinguishes "marginal" from "egregious" poses.

**Goal:** Quantify how far each docked pose is from its lowest-energy free conformer in MMFF94 energy units.

**Approach:** Generate a reference conformer ensemble (ETKDGv3 + MMFF94), make the docked and reference molecules chemically consistent by adding explicit hydrogens to both, relax only the added docked-pose hydrogens while fixing all heavy atoms, take the lowest sampled reference energy as baseline, and report `docked_energy - min_ref_energy` as a relative strain diagnostic. This is not a rigorous solution-phase conformational free energy.

Run `ligand_strain_mmff(docked_sdf, n_ref_conf=20)` from `scripts/validate_poses.py`; it returns a DataFrame with `pose_idx`, `strain_kcal` and a `note` for poses it could not score.

Interpret relative MMFF strain in the context of ligand chemistry, conformer-sampling coverage, and force-field support. Boström et al. (1998) found a conformational energy penalty of no more than 3 kcal/mol for about 70% of 33 protein-bound ligands; that result does not establish a universal acceptance cutoff. Treat unusually high values as a prompt for inspection or use a project-defined threshold validated for the series.

## vdW Overlap with Protein

The 2024 benchmark criterion limits protein-ligand overlap to 7.5% of the ligand vdW volume, using protein radii scaled by 0.8. PoseBusters' `bust(...)` computes this check; do not substitute an unvalidated pairwise-distance sketch for its volume calculation.

## Aromatic Ring Planarity

Run `python scripts/aromatic_planarity.py poses.sdf` (or import `aromatic_planarity(mol)` from it); it prints the maximum out-of-plane deviation per pose.

This reimplementation approximates PoseBusters' internal flatness computation; it does not reproduce it exactly. Verified against installed PoseBusters 0.6.5 (`aromatic_ring_flatness`) on a displacement series of the same pose: this formula's deviation reliably passes at 0.29 Å and reliably fails at 0.45 Å and above -- roughly double the 0.25 Å figure sometimes quoted for this check. Treat `bust()`'s own `aromatic_ring_flatness` column as authoritative for pass/fail; use this snippet only as a supplementary diagnostic, not a gate.
