# Metabolomics Pipeline Usage Guide

## Overview

This workflow orchestrates an untargeted LC-MS metabolomics study from raw mzML to enriched pathways, chaining five component skills: xcms feature extraction, QC/drift/normalization, confidence-stratified annotation, permutation-validated statistics, and background-aware pathway mapping. It is a sequencer with honest handoffs -- each stage defers to its component skill for parameters and traps, and the workflow's job is to keep the failure modes from one stage from silently corrupting the next. Stable-isotope flux tracing is a separate branch (see metabolomics/isotope-tracing), not part of this pipeline.

## Quick Start

Tell your AI agent what you want to do:
- "Run the full untargeted LC-MS metabolomics pipeline on my mzML files"
- "Process my LC-MS data, correct for drift, and find differential metabolites"
- "Take my feature table through QC, statistics, and pathway mapping"
- "I have no compound IDs -- run mummichog pathway analysis on my m/z peaks"

## Example Prompts

### End-to-End
> "I have centroided mzML files plus pooled QCs from an untargeted study; run xcms preprocessing, QC-correct, find differential metabolites, and map pathways."

> "Process my LC-MS data with the modern xcms 4.x API and report differential metabolites with honest annotation confidence."

### QC and Normalization
> "Correct injection-order drift against my QC samples and confirm it did not absorb biological signal."

> "Filter my feature table by QC RSD and D-ratio, PQN-normalize, and impute the residual missing values by mechanism."

### Statistics and Pathways
> "Run a permutation-validated OPLS-DA and a univariate FDR analysis, then reconcile them."

> "My features have no compound IDs -- run mummichog pathway analysis using the full feature table as background."

## When to Use This Pipeline

- Untargeted LC-MS/MS metabolite profiling
- Metabolic biomarker discovery and treatment-response studies
- Lipidomics (adjust peak widths and annotation; see metabolomics/lipidomics)
- Studies where the feature table must reach pathway interpretation honestly

Not for: stable-isotope flux/fluxomics (see metabolomics/isotope-tracing) or absolute targeted quantification as the primary goal (see metabolomics/targeted-analysis).

## Related Skills

- metabolomics/xcms-preprocessing - Stage 1 feature extraction parameters and the feature-table-as-artifact framing
- metabolomics/normalization-qc - Stage 2 drift correction, RSD/D-ratio filtering, PQN, mechanism-aware imputation
- metabolomics/metabolite-annotation - Stage 3 MSI/Schymanski confidence levels
- metabolomics/statistical-analysis - Stage 4 permutation-validated multivariate and dependence-aware FDR
- metabolomics/pathway-mapping - Stage 5 ORA vs mummichog and background construction
- metabolomics/msdial-preprocessing - Alternative front end entering at Stage 2
- metabolomics/lipidomics - Lipid-specific peak widths and annotation
- metabolomics/targeted-analysis - Absolute quantification branch
- metabolomics/isotope-tracing - Separate stable-isotope flux branch, not a stage of this untargeted pipeline
- multi-omics-integration/mofa-integration - Integrating the feature table with other omics layers
