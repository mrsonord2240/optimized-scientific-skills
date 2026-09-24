# XCMS Untargeted LC-MS Preprocessing Usage Guide

## Overview

This skill drives programmatic untargeted LC-MS feature extraction in R with the modern xcms 4.x API: raw centroided mzML in, a features-by-samples intensity table out, with the traps that make untargeted tables untrustworthy (silently deleted peaks, fabricated gap-fill intensities, mis-registered alignment, one-compound-many-features redundancy) handled in `SKILL.md`. Install commands, prerequisites and the workflow are in `SKILL.md` (Version Compatibility, Decision Tree, and the per-step sections).

## Quick Start

Tell your AI agent what you want to do:
- "Process my centroided mzML files with xcms into a feature table"
- "Detect peaks with CentWave using parameters appropriate for my Orbitrap UHPLC data"
- "Align retention times with obiwarp and group peaks across samples"
- "Gap-fill, then report the fraction of filled values per feature"
- "Collapse adduct and isotope redundancy with CAMERA before annotation"
- "Filter features by QC CV and D-ratio with xcms filterFeatures"

## Example Prompts

### Peak Detection
> "Set up CentWaveParam for qTOF UHPLC data and explain why ppm should not be the spec mass accuracy."
> "My trace metabolites are missing from the table - which of prefilter, noise, and snthresh should I lower first?"
> "This is an older low-resolution quadrupole run - what xcms peak picking should I use?"

### Alignment and Grouping
> "Align retention times to a pooled QC reference rather than the first file, then regroup."
> "Choose between obiwarp and peakGroups for a heterogeneous case/control cohort with few shared peaks."
> "Set the grouping bw relative to my residual post-alignment RT scatter on UHPLC."

### Gap-Filling and QC
> "Gap-fill with ChromPeakAreaParam but flag the filled values and report per-feature filled fraction."
> "Filter features with RsdFilter and DratioFilter using my QC and study sample indices."

### Redundancy and Export
> "Run CAMERA in the correct order to collapse adducts and isotopes before annotation."
> "Export the feature table with mz, rt, and intensities for downstream statistics."

## Related Skills

- metabolomics/normalization-qc - Drift correction, CV/D-ratio filtering, and feature-table normalization
- metabolomics/metabolite-annotation - Identification of features into named metabolites
- metabolomics/msdial-preprocessing - GUI/MS-DIAL alternative with MS2Dec deconvolution and GC-EI support
- metabolomics/statistical-analysis - Differential and multivariate statistics on the feature table
