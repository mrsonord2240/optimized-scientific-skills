# Combinatorial CRISPR Screens - Usage Guide

## Overview

Decision-grade design and analysis of combinatorial CRISPR screens. Covers paired-Cas9 (Big Papi, Najm 2018), enhanced AsCas12a multiplex (enCas12a, DeWeirdt 2021), in4mer 4-guide arrays (Esmaeili Anvar 2024), and Inzolia paralog library; paralog buffering detection (Dede 2020, Thompson 2021) for synthetic-lethal targets; genetic-interaction (GI) scoring; and the MAGeCK MLE interaction-term approach.

## Prerequisites

```bash
conda install -c bioconda mageck   # not on PyPI
pip install pandas numpy scipy matplotlib statsmodels   # statsmodels for GI BH-FDR correction
# Inzolia library annotation: bioconductor or hart-lab
# Inzolia annotation: Nat Commun 15:3577 supplement, or Addgene #209551-2

# Cas12a-aware analysis tools
# No dedicated Cas12a helper package; use pandas/numpy with the library annotation
```

Required inputs:
- Cassette-level count matrix (rows = cassettes, columns = samples)
- Library annotation: cassette_id -> gene_A, gene_B (and gene_C, gene_D for in4mer)
- Singleton controls in library
- Vehicle / Day 0 / endpoint samples

## Quick Start

Tell the AI agent what to do:
- "Design an Inzolia paralog-pair library covering 600 paralog pairs in the kinase / transcription factor / DDR families with all-singleton controls"
- "Analyze my Cas12a multiplex paralog screen: compute per-pair GI scores; identify synthetic-lethal pairs at GI z-score <-2"
- "Compare Big Papi paired-Cas9 vs Inzolia Cas12a for the same paralog-pair set"
- "Diagnose: my paired-Cas9 cassettes seem to act as single perturbations -- diagnose dual-cassette cloning failure"
- "Run MAGeCK MLE with explicit interaction terms on my 4-condition combinatorial screen (NT, A_KO, B_KO, A_B_KO)"

## Example Prompts

### Library Design

> "Build an Inzolia-style enAsCas12a 4-guide array library covering 400 receptor tyrosine kinase paralog pairs. Include singletons (gene A alone, gene B alone, double-NTC) so we can compute GI = double_LFC - sum(single_LFC). Output library.tsv."

> "Design Big Papi constructs for 150 specific DDR pathway pairs using the pPapi architecture: SaCas9 sgRNA and SpCas9 sgRNA expressed from the U6 and H1 promoters, with SpCas9 supplied by the cell line. Verify by amplicon sequencing of clones."

### Per-Pair GI Analysis

> "Apply GI scoring to my paired-Cas9 screen. From paired_lfc.tsv and single_lfc.tsv, compute observed_double - (single_A + single_B). Z-normalize across pairs. Classify: SL (z <-2), rescue (z >2), no-interaction (otherwise)."

> "From my in4mer Inzolia screen counts, aggregate cassette-level LFCs by paralog pair. Compute per-pair GI z-scores. Identify 24+ synthetic-lethal pairs as in Dede 2020."

### MAGeCK MLE with Interactions

> "Build MAGeCK MLE design matrix for a 4-condition combinatorial screen: NT, A_KO, B_KO, A_B_KO. Include interaction column. Run mageck mle; output per-pair `interaction|beta` and `|fdr`. Flag synthetic-lethal pairs at FDR <0.05 with negative beta."

### Method Comparison

> "Compare GI z-scoring vs MAGeCK MLE interaction-term approach on the same screen. Investigate disagreement: which method is more sensitive for small-effect interactions?"

> "Validate Big Papi (Cas9) vs Inzolia (Cas12a) on the same 50 paralog pairs. Quantify concordance of synthetic-lethal calls."

### Diagnostics

> "My Cas12a in4mer screen has 20% of cassettes editing <50% per Cas12a editing assay. Investigate: locus-specific inefficiency? Filter and re-analyze."

> "Many GI z-scores are positive (synthetic rescue) for essential pairs (e.g., RPS3 + RPL11). Diagnose: linear-space saturation; switch to log-space (LFC) GI scoring."

> "My paired-Cas9 cassettes show as single perturbation despite design. Diagnose: recombination between repeated U6/tracrRNA elements; re-clone in pPapi (orthologous SaCas9 + SpCas9, U6 + H1); re-clone with validated dual-cassette protocol."

## What the Agent Will Do

1. Identify screen architecture: paired-Cas9 vs Cas12a multiplex vs in4mer
2. Verify library composition: paired cassettes + singletons + NTCs
3. Run standard pooled-screen counting and QC (see [[screen-qc]])
4. For Cas12a: verify per-locus editing efficiency >50%
5. Compute per-cassette LFCs vs control
6. Aggregate to per-pair LFCs across replicate cassettes
7. Extract single-gene LFCs from singleton-control cassettes
8. Compute GI = observed_double - expected_additive
9. Z-normalize GI across all pairs; apply z = -2 cutoff for synthetic lethal
10. Cross-validate via MAGeCK MLE interaction-term analysis
11. Annotate hits against COSMIC, paralog databases, and known synthetic-lethal pairs
12. Output: per-pair GI table, synthetic-lethal candidates, MAGeCK MLE comparison, validation strategy

## Tips

- Architecture choice (Big Papi vs Inzolia/in4mer vs Perturb-seq): see SKILL.md's Combinatorial Architecture Decision Tree.
- Library/cassette/singleton counts, editing-efficiency and GI-score cutoffs: see SKILL.md's Quantitative Thresholds table.
- Common failure modes (recombined dual-sgRNA constructs, low Cas12a editing efficiency, missing singletons, linear-space saturation, library skew): see SKILL.md's Failure Modes section.
- Drug-target validation of synthetic-lethal hits: see SKILL.md's Cross-Modality Validation section (also covers the research-lead-vs-treatment scope boundary).
- For combinatorial Perturb-seq (single-cell readout of multi-perturbation), see [[perturb-seq-analysis]].

## Related Skills

- crispr-screens/library-design - Inzolia / in4mer / Big Papi design
- crispr-screens/screen-qc - Cassette-level QC
- crispr-screens/mageck-analysis - MAGeCK MLE with interaction terms
- crispr-screens/hit-calling - Cross-method validation
- crispr-screens/perturb-seq-analysis - Combinatorial Perturb-seq
- crispr-screens/copy-number-correction - Pre-correction for cancer-line screens
- pathway-analysis/go-enrichment - Functional analysis of GI clusters
