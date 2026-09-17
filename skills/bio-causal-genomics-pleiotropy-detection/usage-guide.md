# Pleiotropy Detection - Usage Guide

## Overview

Detect and adjust for horizontal pleiotropy in two-sample Mendelian randomization by distinguishing uncorrelated (UHP) from correlated (CHP) pleiotropy and choosing a method battery whose assumptions span both regimes. See SKILL.md's "UHP vs CHP: The Central Postdoc-Grade Distinction" for the regime table, the InSIDE assumption, and the operational rg-threshold rule that triggers a CHP-aware method (CAUSE, LHC-MR, or LCV).

## Operational Decision Flow

See SKILL.md's "Operational Decision Flow (4 Steps)" for the full escalation logic (LDSC rg gate, standard battery, CHP escalation, triangulation). Not repeated here so the two files can't drift apart on a fix like the LCV field-name correction above.

## Prerequisites

R packages install from CRAN and GitHub; see SKILL.md's Version Compatibility section for the exact install commands, version pins, and expected input formats.

## Quick Start

Tell your AI agent what you want to do:
- "Run the standard MR sensitivity battery on my TwoSampleMR-harmonized data"
- "I suspect a shared heritable confounder; run CAUSE alongside IVW, Egger, and MR-PRESSO"
- "My instruments are mostly weak (mean F < 20); use MR-RAPS instead of IVW"
- "I have a polygenic exposure with few significant SNPs; run LHC-MR or LCV"
- "Cluster my instruments by causal estimate to identify heterogeneous mechanisms"
- "Run bidirectional MR with Steiger pre-filtering on both directions"
- "Apply SIMEX to my MR-Egger because I^2_GX is 0.75"
- "Build a STROBE-MR sensitivity reporting table from my MR results"

## Example Prompts

### Standard sensitivity battery (UHP-focused)

> "Run IVW, MR-Egger, weighted median, weighted mode, Cochran Q, Egger intercept, MR-PRESSO with 5000 distributions, Steiger directionality, and leave-one-out on my harmonized data and produce a comparison table."

> "Compute I^2_GX for my Egger analysis; if below 0.9 apply SIMEX correction."

### Suspected CHP (correlated pleiotropy)

> "LDSC genetic correlation between my exposure and outcome is 0.45. Run CAUSE in addition to MR-PRESSO and compare delta_ELPD against the sharing model."

> "Use LHC-MR to jointly estimate forward causal effect, reverse causal effect, and heritable-confounder contribution from genome-wide sumstats."

> "Compute LCV gcp from LDSC-merged sumstats to distinguish causation from pure genetic correlation."

### Weak-IV regime

> "Mean F-statistic across my instruments is 14. Run MR-RAPS with Huber robust loss and overdispersion modeling; report point estimate and CI alongside IVW."

> "Apply MR-Mix and contamination mixture to my weak-IV dataset; reconcile against MR-RAPS."

### Heterogeneous mechanisms

> "I suspect my LDL instruments operate through multiple lipoprotein pathways. Run MR-Clust and report per-cluster IVW estimates with biological annotation suggestions."

### Polygenic exposure

> "My exposure is polygenic with only 40 genome-wide-significant SNPs. CAUSE is underpowered. Run LHC-MR or report LCV gcp instead."

### Drug-target cis-MR

> "Run colocalization at the cis-locus instead of Egger because I only have 3 SNPs in tight LD; complement with Steiger directionality."

### Bidirectional MR

> "Instrument both directions and apply Steiger pre-filter before primary IVW; report forward and reverse IVW + Egger + median + mode side by side."

### SIMEX rescue for Egger

> "My I^2_GX is 0.78. Apply SIMEX correction to MR-Egger using the simex package; report SIMEX-corrected slope alongside naive slope."

### STROBE-MR reporting

> "Generate a STROBE-MR compliant table summarizing all methods, F-statistics, I^2_GX, Egger intercept, PRESSO global / distortion / corrected, Steiger, and CAUSE delta_ELPD with citations."

## What the Agent Will Do

See SKILL.md's "Operational Decision Flow (4 Steps)" and "Standard Sensitivity Battery (Working Reference)" for the full step-by-step behavior, from instrument verification through STROBE-MR table output.

## Related Skills

See SKILL.md's Related Skills section.
