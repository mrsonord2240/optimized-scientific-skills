# Peptide Identification - Usage Guide

## Overview
Match MS/MS spectra to peptide sequences by database search (or spectral-library search), then control false discovery rate with target-decoy competition. The deliverable an agent should act on is a q-value (list-level error) or PEP (per-ID error), never a raw engine score, and PSM-level FDR is a separate problem from protein-level FDR. `SKILL.md` holds the install notes, the command-line route, the traps and the thresholds; this guide is only for choosing the Skill and phrasing the request.

## Quick Start
Tell your AI agent what you want to do:
- "Build a concatenated target-decoy FASTA from my UniProt download and check the decoy count"
- "Run a concatenated target-decoy database search on my mzML against UniProt human"
- "Search this mzML with Sage and rescore the pin with Percolator to a 1% PSM list"
- "Configure trypsin with 2 missed cleavages, 10 ppm precursor and 0.02 Da fragment tolerance"
- "Annotate decoys and filter PSMs to 1% FDR with a proper q-value"
- "Explain why my PEP <= 0.01 cutoff kept so many fewer peptides than q <= 0.01"

## Example Prompts

### Database Search Setup
> "Configure a database search with trypsin, 2 missed cleavages, carbamidomethyl C fixed and oxidation M variable"

> "Route an open-search request for unknown modifications to ptm-analysis"

### Running Searches
> "Run a peptide search against a concatenated target-decoy human FASTA with pyOpenMS"

> "Pick a search engine for varied fragmentation across instruments and justify the choice"

### FDR and Rescoring
> "Annotate target/decoy from the DECOY_ prefix and filter to 1% peptide FDR"

> "Compute q-values from this PSM table assuming a concatenated competition search"

> "Rescore the search with Percolator and report the ID gain at q <= 0.01"

### Results Processing
> "Explain the difference between PEP and q-value for these PSMs"

## Related Skills
- protein-inference - Group peptides to protein groups and control protein-level (picked) FDR
- ptm-analysis - Open/variable-mod search follow-up and per-site PTM localization
- dia-analysis - DIA peptide-centric extraction and scoring; entrapment FDR validation
- quantification - FDR-filtered IDs feed label-free/TMT intensity quantification
- spectral-libraries - Empirical and predicted spectral-library search as an ID alternative
- data-import - Load mzML/raw MS data before identification
- database-access/uniprot-access - Build the target FASTA (canonical vs isoform, contaminants)
