# Analytical Validation and Detection Limits - Usage Guide

## Overview

Establishes and audits liquid-biopsy sensitivity claims using explicit molecule
counts, an empirical dilution-series LoD95 with uncertainty, and a strict
boundary between theoretical sampling calculations and achieved assay
performance.

## Prerequisites

```bash
pip install numpy scipy statsmodels
```

## Quick Start

Tell your AI agent what you want to do:

- "How many genome equivalents give 95% probability that a 0.1% variant is physically present?"
- "Estimate LoD95 and its confidence interval from my replicated binary dilution series."
- "Compute the theoretical sampling-only k-of-N bound for a 30-variant MRD panel, and list every assumption."
- "Audit this 'detects 0.1% VAF' claim for missing input mass, recovery, consensus depth, background, and replicate evidence."

## Interpretation Boundaries

- The default conversion is 303.03 haploid GE/ng, derived directly from 3.3 pg per haploid genome. State and pass `ge_per_ng` if the laboratory uses another convention.
- A Poisson result is template-presence probability only. It excludes recovery, consensus depth, background error, false positives, and calling behavior.
- The multi-locus binomial result is an optimistic theoretical sampling VAF lower bound, not panel LoD95.
- An empirical LoD95 fit needs replicated binary outcomes at three or more VAF levels, observed rates bracketing 95%, an identified positive slope, convergence diagnostics, and a confidence interval.
- `simulate_dilution_series` is a contrived probit teaching model whose `true_lod_vaf` sets its 95% point. It is not assay-validation evidence.

## Reference-Material Identity

- HCC1395 tumor (ATCC CRL-2324) and HCC1395BL matched normal (ATCC CRL-2325) underpin the distinct SEQC2 somatic call set v1.2 (SRA SRP162370; Fang 2021, DOI 10.1038/s41587-021-00993-6).
- SEQC2 Sample A is instead an equal-mass pool of ten UHRR cancer-cell-line DNAs. Sample B is Agilent male reference DNA 5190-8848; their diluted, fragmented Df/Ef/Ff materials are the liquid-biopsy proficiency materials (BioProject PRJNA677999; DOI 10.1038/s41587-021-00857-z).
- Choose the material for the validation question and record the exact material, dilution, fragmentation state, truth set, and commutability caveat.

## Related Skills

- ctdna-mutation-detection - applies these limits to low-VAF somatic calls
- longitudinal-monitoring - per-timepoint LoD and left-censoring of undetectable samples
- tumor-fraction-estimation - the ~3% CNA-based detection floor as an LoD
- experimental-design/multiple-testing - repeated-surveillance specificity and FDR
- clinical-biostatistics/power-and-sample-size - validation-study design
