# Metabolite Annotation Usage Guide

## Overview

Metabolite annotation turns untargeted LC-MS/MS features (m/z, RT, MS/MS) into named compounds, but the central job is honesty: every name must carry a confidence level. This skill guards against the field's recurring failures -- reporting a database hit as an identification, treating a high cosine score as proof, claiming a specific isomer MS/MS cannot resolve, and letting ambiguous annotations poison downstream pathway analysis.

## Prerequisites

Install notes and tool versions are in SKILL.md (Version Compatibility). You need a feature table from upstream preprocessing (metabolomics/xcms-preprocessing or metabolomics/msdial-preprocessing) with ion families already collapsed, and the ion mode and expected adducts.

## Quick Start

Tell your AI agent what you want to do:
- "Match my MS/MS spectra against a reference library and report the matched-peak count, not just the score"
- "Run SIRIUS for molecular formula and compound class on features with no library spectrum"
- "Assign a defensible MSI/Schymanski confidence level to each annotation given the evidence I have"
- "Tell me which annotations are safe to carry into pathway analysis and which are Level 3 hypotheses"

## Example Prompts

### Spectral Library Matching
> "Score my query MS/MS against MassBank using modified cosine with a 0.7 score and 6-peak floor."
> "Use spectral entropy similarity for identity matching and flag anything below the 0.75 natural-products threshold."

### In-silico Annotation
> "Run the SIRIUS 6 subcommand chain for formulas, fingerprints, structures against the bio database, and CANOPUS."
> "I trust formula more than structure -- report the ZODIAC formula and only call a structure confident if COSMIC FDR is set."

### Confidence Assignment
> "Given a library match but no in-house standard, what Schymanski level is this and why?"
> "Collapse my evidence set into a single confidence level and explain what would promote it."

### Avoiding Over-claiming
> "Check whether any of these annotations claim a specific isomer that MS/MS cannot resolve."
> "Before pathway analysis, flag features whose ambiguous candidate sets would inflate enrichment."

The method, tool choice, thresholds and failure modes are in SKILL.md.

## Related Skills

- metabolomics/xcms-preprocessing - Upstream feature extraction
- metabolomics/msdial-preprocessing - Alternative feature extraction and deconvolution
- metabolomics/pathway-mapping - Downstream enrichment that must respect these confidence levels
- metabolomics/lipidomics - Lipid-specific annotation and structural resolution
- proteomics/spectral-libraries - Related spectral-matching concepts
