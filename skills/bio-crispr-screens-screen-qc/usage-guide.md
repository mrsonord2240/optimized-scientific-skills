# CRISPR Screen QC - Usage Guide

## Overview

Decision-grade quality control for pooled CRISPR screens. Covers six bottleneck stages (plasmid pool, Day-0 infection, selection, endpoint, biological signal, copy-number artifact), with stage-specific Gini, skew, replicate Pearson/Spearman, sequencing depth, MOI verification, and essentialome PR-AUC against CEGv2 (Hart 2017). Outputs a DepMap-style composite score and recommends remediation for each failing stage.

## Prerequisites

Install commands and required inputs: see SKILL.md, "Install and Inputs".

## Quick Start

Tell the AI agent what to audit:
- "Audit my Brunello-screen counts: plasmid Gini, Day-0 to endpoint dropout, replicate Pearson and Spearman, CEGv2 PR-AUC, depth in reads per sgRNA, MOI verification"
- "Diagnose why hits include ERBB2 in HER2+ SK-BR-3 -- is this copy-number bias or real essentiality"
- "Decide whether my screen passes DepMap quality thresholds, fails, or is salvageable"
- "Recommend which hit-calling method (MAGeCK RRA / MLE / BAGEL2 / Chronos / drugZ) my screen quality grade supports"
- "Generate a composite quality score across all my replicate pairs to gate which conditions enter hit calling"

## Example Prompts

### Library and Plasmid Pool QC

> "Compute Gini coefficient and skew ratio (p90/p10) on my plasmid-pool sequencing. Pass at Gini <0.1 (MAGeCK-VISPR), zero-count <0.5% (Joung 2017), and skew <2 as a stricter modern convention than Joung 2017's <10. Decide whether to proceed."

> "My plasmid Gini is 0.18 and skew is 4.2. Diagnose: PCR over-amplification, synthesis defect, or cloning bottleneck? Recommend remediation."

### Replicate and Depth Audit

> "Compute pairwise Pearson on log10(counts+1) and Spearman on raw ranks between all replicates within each condition. Flag pairs below Pearson 0.85 or Spearman 0.7."

> "One of my treatment replicates shows Pearson 0.78 vs the other two replicates which are >0.95 with each other. Decide whether to drop or rescue the outlier replicate."

> "Verify sequencing depth: reads per sgRNA per sample. Fail below 100 (Joung 2017 plasmid-QC floor), warn between 100 and 300 (MAGeCK-VISPR), pass at 300+, and treat 500+ as ideal for screening (Joung 2017)."

### Biological Signal (Essentialome Recovery)

> "Compute PR-AUC against CEGv2 essentials (Hart 2017) and NEGv1 non-essentials (Hart 2014). My screen passes only if PR-AUC >0.7. Tell me whether this screen has interpretable biology before I run hit calling."

> "PR-AUC is 0.45 even though Gini and Pearson pass. Diagnose: Cas9 not selected pre-screen, premature timepoint, or library targets wrong TSS?"

### Copy-Number Artifact Diagnostic

> "Run the copy-number bias diagnostic. Compute Spearman ρ between gene-level LFC and copy number from the matched WGS profile. Flag CN bias if abs(ρ) >0.1 and recommend CRISPRcleanR / CERES / Chronos correction."

> "Hits include ERBB2 in SK-BR-3 and FGFR1 in head-and-neck lines. Confirm whether these are copy-number artifacts before publication."

### Composite Quality Gate

> "Generate the DepMap-style composite quality score across all metrics: Gini-inverse, replicate Pearson minimum, PR-AUC, log-depth, detected-fraction. Use this to decide which conditions enter downstream hit calling."

### MOI Verification

> "Verify infection MOI: from the titration plate (8 wells with 1:2 serial dilution), interpolate the infection efficiency at the volume used in the screen. Compute Poisson P(≥2 sgRNAs/cell) at the resulting MOI."

## Related Skills

- crispr-screens/library-design - Design library with adequate skew margin
- crispr-screens/mageck-analysis - Generate count files for QC
- crispr-screens/copy-number-correction - Remediate CN artifact
- crispr-screens/batch-correction - Address inter-batch drift
- crispr-screens/hit-calling - Pick analysis method by quality grade
- crispr-screens/in-vivo-screens - In-vivo-specific bottleneck QC
- read-qc/quality-reports - General NGS QC upstream
