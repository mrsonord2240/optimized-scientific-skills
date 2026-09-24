# Lipidomics Usage Guide

## Overview

Lipidomics is a combinatorial-structure problem wearing a quantification problem's clothes. This skill keeps the agent honest about two things software routinely overstates: the structural-resolution level a lipid name actually carries (the shorthand separator `space` -> `_` -> `/` -> `(9Z)` encodes what was measured, and sn-position is almost never measured), and quantification (one isotope-labeled internal standard per class is non-negotiable because ESI response is head-group-dependent). It guards against in-source-fragment phantom lyso-lipids, sn over-claims, ether/plasmalogen ambiguity, and invalid cross-class comparisons.

See SKILL.md for install, inputs and the workflow the agent follows.

## Example Prompts

### Annotation honesty
> "Parse these lipid names through Goslin and report the structural-resolution level each one actually claims."
> "Re-emit my LipidSearch output at molecular-species level since we only ran CID - drop the sn slashes."
> "Flag any plasmalogen (P-) call that lacks vinyl-ether diagnostic evidence."

### Quantification
> "Normalize each lipid class to its matched internal standard, not a single global standard."
> "Is comparing PE to PC abundance valid in my dataset given the standards I used?"
> "Set up class-based quantification using my EquiSPLASH internal standards."

### Differential and enrichment analysis
> "Run differential lipid analysis between treatment and control and make a class-faceted volcano plot."
> "Test whether any lipid class, chain length, or unsaturation pattern is enriched among the changed lipids."

### Artifact triage
> "My LPC pool is unexpectedly high - check retention-time co-elution against the parent PCs."
> "This apparent odd-chain PC 33:1 - is it real or an isotope/in-source artifact?"

## Related Skills

- metabolomics/xcms-preprocessing - Upstream peak detection and feature extraction
- metabolomics/msdial-preprocessing - MS-DIAL alignment and deconvolution upstream of lipid annotation
- metabolomics/metabolite-annotation - General (non-lipid) annotation and confidence levels
- metabolomics/normalization-qc - Sample normalization and QC framing
- metabolomics/statistical-analysis - Multivariate stats on the lipid abundance matrix
