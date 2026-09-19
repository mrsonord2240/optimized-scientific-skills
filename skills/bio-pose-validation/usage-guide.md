# Pose Validation Usage Guide

## Overview

Validate docked or AI-generated protein-ligand poses for physical plausibility using PoseBusters. Quantitatively measure ligand strain, geometric distortion, vdW overlap, and stereochemistry preservation. Report dataset-specific PB-valid rates rather than assuming a fixed failure rate for an AI docking method. See `SKILL.md` for install, the PoseBusters check/config tables, code patterns, and the Common Errors table.

## Example Prompts

### Standard pose QC
> "Run PoseBusters dock config on docked_poses.sdf with receptor.pdb. Output a table per pose with each check pass/fail. Filter to PB-valid; rank by Vina score."

### Strain energy quantification
> "For each docked pose, compute relative MMFF94 strain versus the lowest sampled reference conformer. Report high-strain outliers for inspection without imposing an unvalidated universal cutoff."

### AI docking validation
> "Run DiffDock-L on 100 ligands. For each, run PoseBusters; report PB-valid rate. Compare to GNINA classical docking on same compounds."

### FEP input prep
> "Filter docked poses to those passing PoseBusters. Report relative MMFF strain and apply only the project-defined cutoff validated for this chemical series before FEP setup."

### Chirality/stereo check
> "Identify chirality-inverted poses in my DiffDock output."

## Related Skills

- chemoinformatics/virtual-screening - Source poses from classical docking
- chemoinformatics/ml-docking-rescoring - DiffDock + GNINA + PoseBusters hybrid
- chemoinformatics/molecular-io - SDF parsing
- chemoinformatics/conformer-generation - Reference conformers for strain
- chemoinformatics/free-energy-calculations - PB-valid input for FEP
- chemoinformatics/covalent-design - Covalent pose validation
