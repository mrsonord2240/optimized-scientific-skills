# Conformer Generation Usage Guide

## Overview

Generate 3D conformer ensembles for molecules from 2D structures. Default to ETKDGv3 (Wang et al. 2020, building on Riniker & Landrum 2015) with MMFF94 optimization for drug-like molecules, and use CREST with GFN2-xTB for difficult macrocycles or peptides. Treat GeoMol as a research alternative whose released-model coverage and target chemistry require validation.

## Prerequisites

See SKILL.md's Installation section (includes the Windows/CREST caveat).

## Quick Start

Tell the AI agent what to do:
- "Generate 50 ETKDGv3 conformers for this SMILES with MMFF94 optimization"
- "Build a 3D conformer ensemble suitable for docking input"
- "Sample macrocycle conformers using CREST + GFN2-xTB"
- "Generate Boltzmann-weighted descriptor averages over a 100-conformer ensemble"
- "Prune conformers below 0.5 RMSD and within a 5 kcal/mol energy window"

## Example Prompts

### Drug-like conformer ensemble
> "Generate 20 ETKDGv3 conformers for each compound in library.csv. MMFF94 optimize. Prune at RMSD < 0.5 A and energy window 10 kcal/mol. Output SDF with conformer index in property."

### Macrocycle high-quality sampling
> "Sample 200 conformers of cyclosporine (CAS 59865-13-3) with CREST + GFN2-xTB. Filter to 5 kcal/mol energy window. Output cluster centroids."

### 3D QSAR descriptor input
> "For each compound in active_set.sdf, generate 50 ETKDGv3 conformers and compute Boltzmann-averaged 3D descriptors (asphericity, eccentricity, PMI). Output descriptor table."

### Docking input
> "Generate a single low-energy ETKDGv3 + MMFF94 conformer for each SMILES; write to multi-SDF for downstream Vina docking."

For the step-by-step workflow, decision tree, n_conf heuristic, and per-method tips, see SKILL.md
(its Decision Tree, Macrocycle Handling, and Common Errors sections cover all of this once).

## Related Skills

- chemoinformatics/molecular-io - Parse molecules
- chemoinformatics/molecular-descriptors - 3D descriptors from ensembles
- chemoinformatics/shape-similarity - Multi-conformer 3D shape searching
- chemoinformatics/virtual-screening - Generate 3D for docking
- chemoinformatics/free-energy-calculations - Conformers for MD/FEP setup
- chemoinformatics/pharmacophore-modeling - 3D pharmacophore from ensembles
