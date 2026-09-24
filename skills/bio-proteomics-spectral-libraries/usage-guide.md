# Spectral Libraries - Usage Guide

## Overview
Builds and manages DIA spectral libraries, which are tables of peptide query parameters (precursor m/z, a few fragment m/z plus relative intensities, normalized RT, optional CCS) rather than whole spectra. Covers experimental DDA, chromatogram, predicted (Koina-served Prosit/AlphaPeptDeep/MS2PIP, DeepLC), and empirically-corrected libraries, plus iRT/CiRT RT calibration, NCE tuning, format conversion, and QC/merge. The point that decides library quality: a predicted library is only as good as its empirical RT/CCS calibration, and NCE must match the predictor's training.

## Prerequisites
Install notes and version checks are in SKILL.md ("Version Compatibility"); the first ms2pip call
downloads a large model, see its Common Errors table.

## Quick Start
Tell your AI agent what you want to do:
- "Generate a predicted library with Prosit via Koina for my protein list"
- "Calibrate my predicted iRT to the observed retention times in my run"
- "Convert my Spectronaut library to DIA-NN format"
- "Merge two libraries without dropping distinct charge states"
- "Build an empirically-corrected predicted library for my non-model organism"

## Example Prompts

### Generating Predicted Libraries
> "Use Prosit via Koina to predict fragment intensities and iRT for this peptide list"

> "Predict a full library including CCS with AlphaPeptDeep for my timsTOF data"

> "Scan candidate NCE values and pick the one that best matches my real spectra"

### Calibrating Retention Time
> "Fit my predicted iRT to observed RT using the Biognosys iRT peptides"

> "I have no iRT spike-in -- calibrate retention time using CiRT endogenous peptides"

> "My gradient is nonlinear; use a LOWESS fit instead of a linear iRT alignment"

### Empirical and Chromatogram Libraries
> "Build a chromatogram library from my GPF-DIA runs with EncyclopeDIA"

> "Build an empirically-corrected predicted library by searching one GPF-DIA pass"

> "Build a DDA library from my FragPipe search results with EasyPQP"

### Format Conversion
> "Convert my library to DIA-NN tsv format and check the RT units"

> "Export my library for OpenSWATH and generate decoys with OpenSwathDecoyGenerator"

> "Reconcile the modification notation between UniMod accessions and delta masses"

### Quality Assessment and Merging
> "Report precursors, proteins, and transitions per precursor in my library"

> "Merge these libraries keeping the full transition key so I don't drop charge states"

## Related Skills

- dia-analysis - Run the DIA search against the library and choose the q-value context
- peptide-identification - Generate the DDA identifications a DDA library is built from
- ptm-analysis - Design PTM-resolved and modified-peptide libraries
- quantification - Summarize and roll up the search output to protein abundances
