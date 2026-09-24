# Quantification - Usage Guide

## Overview
Reconstruct a protein-by-sample abundance matrix from mass spectrometry signals, choosing a summarizer and normalizer that match where the signal physically came from. Label-free (LFQ/MaxLFQ), isobaric (TMT/iTRAQ reporter ions), and metabolic (SILAC) approaches each carry an irreducible error set by their measurement physics, and the peptide-to-protein summarization choice changes the answer more than the downstream statistical test does.

## Quick Start
Tell your AI agent what you want to do:
- "Summarize my MaxQuant peptides to protein level with MSstats"
- "Run the real MaxLFQ algorithm on my peptide intensity matrix"
- "Extract TMT reporter ions and correct for isotope impurity"
- "Bridge my multiple TMT plexes with an IRS reference channel"
- "Compute SILAC heavy/light ratios and check for Arg-to-Pro conversion"

## Example Prompts

### Label-Free Summarization
> "Convert my MaxQuant evidence.txt into normalized protein-level abundances using MSstats Tukey median polish"

> "Run iq::maxLFQ on my peptide quant matrix instead of median centering"

> "Compare TMP and MaxLFQ summarization on the same data and report where the answer moves"

### Normalization
> "Median-center my label-free intensity matrix to correct sample loading"

> "Apply sample-loading normalization then IRS to bridge my three TMT plexes"

### TMT/iTRAQ Processing
> "Extract TMT10 reporter ions from my mzML and apply lot-specific impurity correction"

> "Explain why MS2 reporter quant compresses my fold changes and whether SPS-MS3 helps"

### SILAC
> "Compute SILAC log2 ratios but keep proteins present only in the heavy channel"

> "Check my SILAC labeling efficiency and flag Arg-to-Pro conversion before trusting ratios"

### AP-MS
> "Score my bait pulldown against the GFP control IPs and tell me which prey are real interactors"

## Related Skills
See the Related Skills section of `SKILL.md`.
