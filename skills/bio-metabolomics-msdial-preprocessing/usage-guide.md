# MS-DIAL Preprocessing Usage Guide

## Overview

MS-DIAL turns LC-MS or GC-MS raw data into an aligned feature table, with its distinctive strength being MS2Dec spectral deconvolution that reconstructs clean MS/MS from chimeric DDA/DIA spectra. This skill covers running MS-DIAL (GUI or the cross-platform console), choosing DDA vs DIA, and importing the alignment-result table into R or Python with honest filtering. The trap it guards against: preprocessing is a modeling step, not a neutral clean-up - the same raw files through MS-DIAL vs XCMS give different feature tables and different marker lists, gap-filled low-Fill% features fabricate intensity, and an annotation tag is not an identification.

## Prerequisites

```bash
# MS-DIAL: https://systemsomicslab.github.io/compms/msdial/main.html (GUI is Windows-only)
# Console (cross-platform): a separately-packaged release asset, MSDIALCUI.exe
#   (confirm the current subcommand list by running it with no arguments -- current
#   5.x builds expose one `lcms` subcommand for both DDA and DIA, plus `gcms`)
# Vendor raw -> mzML (DDA/GC) via ProteoWizard msconvert; -> ABF (DIA) via the Reifycs ABF converter
pip install pandas numpy   # for parsing the export in Python
```

Conceptual prerequisites: the acquisition mode (DDA vs DIA/SWATH, set per file rather than by console subcommand) and the platform (LC-ESI vs GC-EI, both supported directly by the current MS-DIAL 5.x console). The MS-DIAL vs XCMS choice should be made before processing.

## Quick Start

Tell your AI agent what you want to do:
- "Run MS-DIAL on my LC-MS DDA folder with the console and produce an alignment table"
- "Process my SWATH/DIA data with MS-DIAL and deconvolve the MS/MS"
- "Parse my MS-DIAL alignment export into a clean feature matrix in R"
- "Filter my MS-DIAL table by Fill% and MS/MS support, mapping annotation tags to MSI levels"
- "Should I use MS-DIAL or XCMS for this cohort?"

## Example Prompts

### Choosing and Configuring
> "I have wide-isolation-window MS/MS data - which MS-DIAL console mode and input format do I need?"
> "Help me set the minimum peak height for an Orbitrap run; the default looks TOF-tuned."
> "I have GC-EI data - can I use MS-DIAL 5, and what alignment scale should I use?"

### Running
> "Build the MSDIALCUI console command to process my DDA folder headless on a Linux cluster."
> "Set the alignment reference to a pooled QC instead of the first file."

### Importing and Filtering
> "Read my MS-DIAL AlignResult export into pandas, skipping the header rows, and split metadata from sample columns."
> "Drop features below 70% Fill and require MS/MS assignment before I trust any identity."
> "Map MS-DIAL annotation tags to MSI confidence levels so I don't overtrust mass-only hits."

## What the Agent Will Do

1. Decide MS-DIAL vs XCMS, and DDA vs DIA vs GC, for the data at hand
2. Build the correct MSDIALCUI console command (`lcms` for DDA/DIA, `gcms` for GC-MS) with the right input format
3. For DIA, or a mixed-mode batch, set `acquisition_type` per file via a CSV passed to `-i`
4. Parse the alignment export (correct header offset, metadata vs sample columns anchored on the header block, not a fixed column position)
5. Filter on Fill% (0-1 fraction), MS/MS support (case-insensitive), and QC/blank thresholds
6. Tie annotation tags to MSI levels and hand off identification to metabolite-annotation

## Tips

- DIA/SWATH mode accepts ABF input only; convert before running, and set `acquisition_type=DIA` per file (it is not a separate console subcommand).
- The GUI is Windows-only and hangs on `--help`; use the separately-packaged MSDIALCUI console for headless and cluster runs.
- Set the alignment reference to a pooled QC, never to file #1 by default.
- A low Fill% means the value is mostly gap-filled noise; report the filled fraction of every hit. MS-DIAL 5.x reports Fill% as 0-1, not 0-100.
- An annotation name without MS/MS is at best a putative (MSI Level 3) ID - require MS/MS before believing it (real value is `True`/`False` text; compare case-insensitively).
- Replicate strong findings across a second pipeline (e.g. XCMS); one-software hits are candidates, not results.

## Related Skills

- metabolomics/xcms-preprocessing - Programmatic R preprocessing and the feature-table-as-artifact framing
- metabolomics/lipidomics - Lipid annotation mode and LipidBlast workflows
- metabolomics/metabolite-annotation - MSI confidence levels and orthogonal-evidence identification
- metabolomics/normalization-qc - Drift correction, QC/CV/D-ratio filtering, MNAR-aware imputation
