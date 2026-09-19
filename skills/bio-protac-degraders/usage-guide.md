# PROTAC and Bivalent Degrader Design Usage Guide

## Overview

Design PROTACs (bivalent molecules recruiting E3 ligase to target for proteasomal degradation). Balance target-ligand, E3-ligand, linker geometry, cooperativity, and cell permeability. Cover ternary complex prediction (PRosettaC, DeepTernary, AlphaFold3), cooperativity, hook effect, and DC50/Dmax.

## Prerequisites

```bash
pip install rdkit numpy scipy
```

`examples/protac_enumerate.py` and `examples/cooperativity_dc50.py` run with the above alone.
`examples/ternary_geometry_screen.py` needs only RDKit too, and is a local linker-reach
pre-filter, not a ternary complex prediction. PRosettaC, DeepTernary, AlphaFold3, Boltz, and
HADDOCK are external services/gated local tools -- see SKILL.md's "Ternary Complex Prediction
Tools" for how to invoke each one; this Skill does not install or run them.

## Quick Start

Tell the AI agent what to do:
- "Design PROTAC for kinase X target using CRBN E3"
- "Enumerate linkers between target-ligand and VHL ligand"
- "Predict ternary complex for kinase + PROTAC + CRBN"
- "Optimize linker length to maximize cooperativity"

## Example Prompts

### CRBN PROTAC design
> "For target kinase X (PDB 5XYZ), design an exploratory PROTAC series using a justified E3 recruiter. Vary linker composition and geometry around candidates compatible with the binary structures. Output connected, synthesis-review-ready SMILES with structural scores reported as hypotheses."

### Ternary complex prediction
> "Predict the ternary structure for this target-ligand-E3 PROTAC complex using PRosettaC. Report structural and interface scores; state that cooperativity alpha requires a binding experiment."

### Linker optimization
> "Given target-ligand and E3-ligand exit vectors in a shared structural hypothesis, suggest a small linker series spanning rigid and flexible chemistries. Explain how conformer feasibility will be checked rather than inferring atom count from distance alone."

### VHL alternative design
> "Switch from CRBN to VHL E3. Adjust linker to maintain ternary geometry. Re-predict ternary complex."

## What the Agent Will Do

1. Define target ligand (from co-crystal or docked) and E3 ligand (pomalidomide / VHL ligand).
2. Compute distance between attachment points on each ligand.
3. Enumerate a linker series around geometries supported by the binary structures.
4. Combine target-linker-E3 SMILES.
5. Predict ternary complex via an external submission (PRosettaC, AlphaFold3, Boltz, or DeepTernary -- see SKILL.md's "Ternary Complex Prediction Tools" for how to invoke each).
6. Score by linker geometry and structural/interface metrics; report experimental cooperativity separately when available.

## Related Skills

- chemoinformatics/molecular-io - Parse ligand SMILES
- chemoinformatics/reaction-enumeration - Linker enumeration
- chemoinformatics/generative-design - REINVENT linker mode
- chemoinformatics/conformer-generation - Ternary conformer sampling
- chemoinformatics/virtual-screening - Validate target ligand binding
- chemoinformatics/free-energy-calculations - Ternary ABFE
- chemoinformatics/admet-prediction - PROTAC ADMET specifics
- structural-biology/structure-io - PDB / mmCIF for ternary complex
