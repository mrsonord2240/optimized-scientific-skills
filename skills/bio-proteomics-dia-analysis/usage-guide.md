# DIA Analysis - Usage Guide

## Overview
Identifies and quantifies proteins from data-independent acquisition (DIA) mass spectrometry by scoring reconstructed fragment-chromatogram peak groups against a decoy null, then filtering the output at the correct q-value level and context. Covers DIA-NN library-free search (which in DIA-NN is a predicted-library search) and library-based routes, with notes on Spectronaut, OpenSWATH, and EncyclopeDIA. The crux: every wide-window MS2 is chimeric, so the engine deconvolves rather than matches, and "1% FDR" only means something once the level (precursor/peptide/protein-group) and context (run vs experiment-wide/global) are named.

## Prerequisites
```bash
pip install pandas pyarrow numpy
# CLI: DIA-NN (recommended), MSFragger-DIA/FragPipe, OpenSWATH, EncyclopeDIA
# Commercial: Spectronaut
# Staggered/overlapping data: DIA-NN supports them natively; other engines need ProteoWizard msconvert (peak picking first, then demultiplexing)
```

## Quick Start
Tell your AI agent what you want to do:
- "Run DIA-NN on my mzML files via the predicted-library route at 1% FDR"
- "Search my DIA data against a spectral library I built"
- "Load the DIA-NN report.parquet, filter q-values correctly, and build a protein matrix"
- "My cohort hits do not reproduce -- which q-value should I be filtering on?"
- "Demultiplex my staggered-window data before searching"

## Example Prompts

### Library-Free and Predicted-Library Analysis
> "Run DIA-NN via the predicted-library route against UniProt human, two-pass, with MS1/MS2 mass accuracy fixed for our timsTOF"

> "Set up library-free directDIA for my Astral 2-Th narrow-window runs"

> "Generate an in-silico predicted library and search my DIA data against it"

### Library-Based Analysis
> "Search my DIA data against the chromatogram library I built in EncyclopeDIA"

> "Run DIA-NN with my Prosit-predicted library for targeted extraction"

> "Use a previously generated .speclib so I do not redigest the FASTA"

### FDR and Output Filtering
> "Filter my DIA-NN report at precursor and protein-group level, and global PG for the matrix"

> "My pg_matrix and globally filtered report have different protein-group counts -- check the DIA-NN log and explain the version-specific filters"

> "My per-run 1% FDR cohort has irreproducible hits -- fix the FDR context"

### Results Processing
> "Load report.parquet, drop unquantified zeros, and log2-transform a MaxLFQ protein matrix"

> "Pivot the filtered long report into a wide protein-by-run matrix"

> "Set up auditable run/experiment/global q-values with OpenSWATH and PyProphet"

## What the Agent Will Do
1. Identify the acquisition design (fixed/variable/staggered windows, diaPASEF, narrow-window Astral) and demultiplex staggered data at conversion if the engine lacks native overlapping-window support
2. Choose the route -- predicted-library / library-free (default) or library-based -- per the decision tree
3. Generate the predicted library from FASTA without raw files, then run DIA-NN (or the chosen engine) against that library with fixed instrument-appropriate mass accuracy and two-pass global FDR
4. Read report.parquet and filter at the correct LEVEL (precursor + protein-group) and CONTEXT (run + global for matrices)
5. Convert unquantified zeros to NA, then hand the log2 matrix to quantification and differential-abundance

## Tips
- The predicted-library route is the modern default; raw directDIA on wide-window data needs careful two-pass FDR.
- "1% FDR" alone is ambiguous; always state the level (precursor/peptide/protein-group) and context (run/experiment/global).
- For cohorts, filter on Global.PG.Q.Value, not the per-run Q.Value, or the experiment-wide error inflates.
- DIA-NN 1.9+ writes report.parquet by default; loaders assuming report.tsv silently break.
- Matrix and report counts are not interchangeable. DIA-NN 2.6.1 can make the matrix larger than a globally filtered report; inspect `report.log.txt` for the version and actual applied filters.
- DIA-NN 2.x requires predicted-library generation and raw-data analysis as separate stages. Do not combine raw files with `--fasta-search --predictor` in one command.
- An EasyPQP `library.tsv` may lack `FragmentCharge`, `FragmentType`, and `FragmentSeriesNumber` required by DIA-NN; prefer a DIA-NN-format export or validate a converter before cohort-scale search.
- Fix mass accuracies per instrument (timsTOF 15/15, Orbitrap Astral 10/4, TripleTOF 20/20 ppm MS2/MS1); --mass-acc 0 is optimised on the first run and reused, so results depend on run order.
- Aim for >= 6 MS2 points across each LC peak; if quant is noisy, the window/cycle design may be the cause.
- Convert DIA-NN's 0 (not-quantified) to NA before log2 or normalization.

## Related Skills

- spectral-libraries - Build experimental, chromatogram, or predicted libraries to search against
- quantification - Normalization, MaxLFQ roll-up, and matrix summarization after filtering
- differential-abundance - Moderated statistical testing of the protein matrix
- proteomics-qc - Per-run ID counts, missing-value rates, and acquisition QC
- data-import - Convert and load raw vendor/mzML data before the search
- peptide-identification - DDA spectrum-to-peptide matching (the non-DIA counterpart)
- workflows/proteomics-pipeline - End-to-end DIA-to-differential-abundance orchestration
