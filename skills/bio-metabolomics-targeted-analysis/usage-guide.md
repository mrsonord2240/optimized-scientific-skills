# Targeted Metabolomics Analysis Usage Guide

## Overview

Targeted metabolomics quantifies a closed, pre-defined panel of known metabolites and reports absolute concentrations with units, using MRM/SRM on a triple quadrupole or PRM on a high-resolution instrument. The central enemy is the matrix effect: co-eluting matrix suppresses analyte ionization, and only a co-eluting stable-isotope-labeled internal standard truly corrects it. This guide covers building a calibration assay, choosing the internal-standard and weighting strategy, confirming identity by ion ratio, and validating to a depth that matches the decision the number supports.

## Prerequisites

```bash
# Skyline (free, vendor-neutral) for transition lists, integration, export
# Download from: https://skyline.ms/
# R for post-export curve fitting and validation metrics
Rscript -e 'install.packages(c("ggplot2", "dplyr"))'
# Optional Python path
pip install pandas numpy scipy
```

Conceptual prerequisites: an authentic reference standard for each analyte (its purity scales every reported number), ideally one stable-isotope-labeled internal standard per analyte, a defined transition list (quantifier + qualifier per analyte with collision energies), and a decision about how much validation the application demands.

## Quick Start

Tell your AI agent what you want to do:
- "Build a 1/x^2-weighted calibration curve and accept it by back-calculated %RE, not R-squared"
- "Normalize each analyte to its co-eluting internal standard and quantify the samples"
- "Confirm identity with the quantifier/qualifier ion ratio and flag interferences"
- "Compute LOD and LLOQ and mark samples below the validated range as not reportable"
- "Help me decide an internal-standard and validation strategy for a clinical assay"

## Example Prompts

### Calibration and Weighting
> "Fit unweighted, 1/x, and 1/x^2 calibration curves and pick the weighting that minimizes low-end back-calculated relative error."
> "Show the per-level %RE for each calibrator and tell me whether the curve passes ICH M10."
> "Set the LLOQ to the lowest calibrator within +/-20% accuracy."

### Quantification
> "Normalize analyte areas to the SIL-IS and back-calculate concentrations from the response-ratio curve."
> "Apply the dilution factor and report concentrations in micromolar, flagging anything below the LLOQ."

### Identity and Quality
> "Compute the qualifier/quantifier ion ratio per sample and flag any outside +/-30% of the calibrator ratio."
> "Check carryover in the blank injected after the top calibrator."
> "Estimate the matrix factor and the IS-normalized matrix factor across matrix lots."

### Strategy and Validation
> "My panel has 40 chemically diverse metabolites and one global IS -- where is my accuracy at risk?"
> "Lay out the ICH M10 parameters I need for a regulated PK assay versus an exploratory study."
> "Should I use a deuterated or 13C internal standard, and what do I verify before trusting it?"

The full workflow (panel/IS/validation-tier selection, calibration, LOD/LLOQ, IS normalization,
ion-ratio confirmation, precision, and report format) and every failure-mode tip live in
SKILL.md -- it is the single source for what the agent does and why.

## Related Skills

- metabolomics/xcms-preprocessing - Upstream feature detection for untargeted discovery before targeted validation
- metabolomics/statistical-analysis - Group comparison and multivariate analysis of quantified concentrations
- metabolomics/isotope-tracing - Stable-isotope tracing and flux (MID), the adjacent discipline this skill hands off to
- metabolomics/normalization-qc - QC-sample-driven drift correction and RSD filtering
- clinical-biostatistics/cdisc-data-handling - Regulated-trial bioanalysis data handling when targeted numbers feed a clinical study
