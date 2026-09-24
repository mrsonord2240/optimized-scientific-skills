# Covalent Inhibitor Design Usage Guide

## Overview

Design covalent inhibitors targeting Cys, Lys, Ser, Thr, Tyr, or Asp residues. Balance warhead reactivity, GSH stability, geometric accessibility, and irreversible vs reversible covalent. Covers DOCKovalent, HCovDock, GOLD covalent, and reactivity-aware SAR.

## Prerequisites

Versions and installs: see SKILL.md's Version Compatibility.

## Quick Start

Tell the AI agent what to do:
- "Suggest warheads for cysteine-selective covalent inhibitor"
- "Score acrylamide-containing analogs for reactivity"
- "Plan reactive group SAR for KRAS G12C inhibitor"
- "Classify candidate warheads, then identify which compounds require experimental GSH-reactivity testing"

## Example Prompts

### Cysteine targeting
> "Identify acrylamide and chloroacetamide candidates in library.smi. Report alpha substitution as a structural feature and flag compounds for matched GSH-reactivity measurements."

### Reactivity SAR
> "For 30 acrylamide compounds in series.csv, compute alpha-C substituent count and compare it with measured GSH rates and kinact/Ki without treating it as a LUMO calculation."

### Covalent docking
> "Dock acrylamide candidates against EGFR C797 with a reaction-appropriate protocol. Inspect Cys797 Sγ-to-electrophile distance and approach geometry, then integrate measured reactivity."

### Bivalent / PROTAC
> "Combine cysteine-warhead inhibitor with E3 ligand via 12-atom linker. Predict ternary complex stability."

## What the Agent Will Do

1. Identify warheads in input (acrylamide, chloroacetamide, vinyl sulfone, etc.).
2. Compute structural features and integrate measured or validated reactivity data.
3. Run covalent docking if requested (DOCKovalent / HCovDock / GOLD).
4. Validate reaction-atom geometry: for cysteine, use Sγ and the ligand electrophilic atom, not Cβ.
5. Review experimental GSH reactivity and off-target evidence.
6. Output ranked candidates with covalent pose + reactivity tier.

## Tips

For reactivity, GSH-stability, geometry, and kinact/Ki caveats, see SKILL.md's Reactive Residue Taxonomy, Intrinsic Reactivity Assays, Reactivity Surrogates, Kinetics: kinact/Ki, and Per-Tool Failure Modes sections — do not restate them here.

## Related Skills

- chemoinformatics/substructure-search - Warhead SMARTS detection
- chemoinformatics/virtual-screening - Non-covalent docking first
- chemoinformatics/pose-validation - Validate covalent docking
- chemoinformatics/molecular-descriptors - Reactivity surrogates
- chemoinformatics/admet-prediction - ADMET of covalent leads
- chemoinformatics/protac-degraders - PROTAC with covalent warhead
