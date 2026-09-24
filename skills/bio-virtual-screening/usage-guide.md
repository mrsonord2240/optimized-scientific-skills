# Virtual Screening - Usage Guide

## Overview
Screen compound libraries against protein targets using molecular docking with AutoDock Vina. Predict binding poses and affinities for hit identification.

## Prerequisites
```bash
pip install vina
pip install rdkit meeko
```

## Quick Start
Tell your AI agent what you want to do:
- "Dock my compound library against this protein"
- "Set up a virtual screen against the binding site"
- "Rank compounds by predicted binding affinity"
- "Prepare my protein for docking"

## Example Prompts

### Receptor Preparation
> "Prepare my protein PDB for docking by reviewing waters, cofactors, metals, and alternate locations; retain justified structural components and add hydrogens."

> "Convert my protein to PDBQT format for Vina."

### Docking
> "Dock this ligand to my protein and show the best poses."

> "Screen my library of 1000 compounds against the active site."

### Analysis
> "Rank my docked compounds by binding affinity."

> "Extract the top 10 hits with affinities better than -8 kcal/mol."

## What the Agent Will Do
1. Prepare receptor (review waters, cofactors, metals, and alternate locations; retain justified structural components; add H and convert to PDBQT)
2. Prepare ligands (generate 3D, minimize, convert to PDBQT)
3. Define binding site from ligand or coordinates
4. Run Vina docking for each compound
5. Collect and rank results by affinity

## Tips
See SKILL.md's "Exhaustiveness", "Box too small" (`references/failure-modes.md`), "Receptor Preparation", "Version Compatibility", and "Common Errors" sections for box sizing, exhaustiveness, receptor-prep, and version-drift guidance -- kept in one place there.

## Related Skills
- molecular-io - Load and convert molecules
- admet-prediction - Filter before docking
- structural-biology/structure-io - Protein structure handling
- structural-biology/modern-structure-prediction - Generate targets
