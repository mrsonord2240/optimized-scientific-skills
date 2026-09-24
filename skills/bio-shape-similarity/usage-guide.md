# Shape Similarity Usage Guide

## Overview

3D shape-based similarity searching using USRCAT (ultrafast), Open3DAlign (RDKit), ROCS (commercial), or ShaEP. Find scaffold-hopped compounds that 2D fingerprints miss; identify bioisosteric replacements via shape + color (Tanimoto-Combo).

## Prerequisites

Install notes are in `SKILL.md` (Version Compatibility).

## Quick Start

Tell the AI agent what to do:
- "Find compounds with similar 3D shape to my query molecule"
- "Score library with USRCAT shape descriptors; top 100 hits"
- "Open3DAlign rescoring on top USRCAT hits"
- "Find scaffold-hopped compounds using shape-high and ECFP4-low cutoffs calibrated on my reference set"

## Example Prompts

### USRCAT pre-filter
> "Compute USRCAT descriptors for query.sdf and library.sdf. Rank library by similarity; return top 500."

### Open3DAlign rescore
> "Take top 500 USRCAT hits and rescore with Open3DAlign for accurate shape alignment. Return top 50."

### Scaffold-hopping
> "Calibrate shape-similarity and ECFP4-dissimilarity cutoffs on my reference set, then output 20 candidate scaffold hops."

### Conformer-aware shape search
> "For each library compound, generate 20 conformers; find best-shape conformer match to query. Use Open3DAlign."

The workflow, thresholds and failure modes are in `SKILL.md`.

## Related Skills

- chemoinformatics/molecular-io - Parse molecules
- chemoinformatics/conformer-generation - Generate conformer ensembles
- chemoinformatics/similarity-searching - 2D similarity comparison
- chemoinformatics/pharmacophore-modeling - Pharmacophore alternative
- chemoinformatics/virtual-screening - Shape as pre-filter
