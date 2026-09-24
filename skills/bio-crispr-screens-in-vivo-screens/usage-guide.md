# In Vivo CRISPR Screens - Usage Guide

## Overview

Decision-grade design and analysis of in vivo CRISPR screens. Covers the bottleneck math (cell number limits force focused libraries); Manguso 2017 immune-evasion methodology; Chen 2015 tumor screens; CRISPR-StAR stochastic post-engraftment sgRNA activation (Uijttewaal 2025 Nat Biotechnol 43:1848) for escaping early bottlenecks; syngeneic vs xenograft vs PDX choice; per-animal clonal variability; tumor-explant DNA extraction; per-animal meta-analysis; in vivo CEGv2 calibration limitations.

## Prerequisites

**Requires prior IACUC (or equivalent institutional) approval before any animal work** — see SKILL.md's "Ethical & Regulatory Requirements" section; this Skill does not itself provide ethical review.

Install notes and versions: SKILL.md "Version Compatibility".

Required inputs:
- Animal cohorts per condition (size: SKILL.md "Quantitative Thresholds")
- Focused library sized to the bottleneck (SKILL.md "Focused Library Design for In Vivo")
- Cas9+ cell line (selected by FACS before infection)
- Plasmid pool sequencing as baseline
- Per-animal tumor DNA + sgRNA amplification primers

## Quick Start

Tell the AI agent what to do:
- "Design a focused in vivo CRISPR library targeting immune-evasion biology following Manguso 2017 methodology: 2,000 kinases + surface proteins + immune factors, 4 sgRNAs/gene"
- "Compute bottleneck math: my B16 syngeneic model can take 2M cells; pick library size for 100x coverage maintainable through bottleneck"
- "Apply CRISPR-StAR so a genome-scale library survives the in vivo bottleneck: activate sgRNAs after engraftment for paired internal controls"
- "Run MAGeCK MLE on my in vivo screen with animal-as-batch covariate; output per-condition beta scores after batch adjustment"
- "Diagnose: why is my in vivo CEGv2 PR-AUC only 0.45 despite passing all in vitro QC?"

## Example Prompts

### Library Design for In Vivo

> "Design a focused library targeting immune-evasion biology in syngeneic B16-OVA mouse melanoma. Include: kinases, cell surface proteins, immune-regulatory genes (Manguso 2017 selection criteria). Total ~2,000 genes x 4 sgRNAs = 8,000 sgRNAs. Verify coverage with 2M cells implanted (= 250x cells/sgRNA before the engraftment bottleneck)."

> "For a syngeneic colorectal cancer model with maximum 3M cells implantable: design library at 5,000 sgRNAs (= 600x effective coverage). Pick genes by pathway relevance (metabolism, immune, proliferation)."

### CRISPR-StAR

> "Set up CRISPR-StAR (Uijttewaal 2025) screen with tamoxifen-inducible CreERT2 sgRNA activation in a syngeneic model. Cells infected with library at MOI 0.3; implanted in 5 mice per condition; tamoxifen induction Day 5 post-implant activates the sgRNA in ~half of each clone for paired internal controls; tumor harvest Day 21. Compute expected coverage maintenance vs a conventional in vivo screen."

### Per-Animal Analysis

> "Run MAGeCK on each animal vs plasmid pool separately; meta-analyze with Stouffer's Z; identify genes consistent across animals. Compare to single combined-animal analysis."

> "Run MAGeCK MLE with animal as batch covariate. Output per-condition beta after batch adjustment. Compare to per-animal RRA + meta-analysis."

### Diagnostics

> "In vivo CEGv2 PR-AUC is 0.45 despite in vitro on same line showing 0.85. Diagnose: context-specific essentialome, clonal dominance, or Cas9 selection failure?"

> "My screen has 100x effective coverage at endpoint but only 60% of library detected. Investigate: bottleneck at engraftment or PCR amplification?"

### Cross-Validation

> "Validate top in vivo hits in vitro (matched cell line, no animal context). Identify hits that are in-vivo-specific (require tumor microenvironment) vs cell-intrinsic."

> "Arrayed validation of top 10 in vivo hits with n=10 mice each. Confirm consistent effect across biological replicates."

## What the Agent Will Do

Follows SKILL.md end to end: model choice, bottleneck-adjusted library size, cohort plan, CRISPR-StAR when genome-scale is needed, tumor DNA recovery, MAGeCK MLE or per-animal RRA + meta-analysis, and the Validation Checklist.

## Related Skills

- crispr-screens/library-design - Focused library design
- crispr-screens/mageck-analysis - MAGeCK MLE with batch covariate
- crispr-screens/hit-calling - Per-animal meta-analysis
- crispr-screens/screen-qc - In vivo-specific QC
- crispr-screens/batch-correction - Animal cohort as batch
- crispr-screens/combinatorial-screens - In vivo combinatorial
- crispr-screens/copy-number-correction - Cancer-line in vivo
- pathway-analysis/go-enrichment - Functional analysis
