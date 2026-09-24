# GSEA - Usage Guide

## Overview
Gene Set Enrichment Analysis (GSEA) tests whether a gene set shifts coordinately toward one end of a genome-wide ranked gene list, without ever choosing a significance cutoff. This Skill covers preranked GSEA with clusterProfiler's gseGO/gseKEGG/gsePathway/GSEA, the matrix-based CAMERA/fry tests that are the correlation-honest default when the expression matrix is available, and per-sample ssGSEA/GSVA scores. Everything the agent acts on (ranking statistic, permutation type, thresholds, failure modes, install) is in `SKILL.md`; a preselected unranked gene list belongs to go-enrichment (ORA).

## Quick Start
Tell your AI agent what you want to do:
- "Run GSEA on my DESeq2 results ranked by the Wald statistic"
- "Find GO biological processes with coordinated expression changes across all my genes"
- "Run GSEA against the MSigDB Hallmark collection on my ranked list"
- "Score each sample for pathway activity with GSVA for downstream clustering"
- "I have the expression matrix and groups: which KEGG pathways differ between case and control?"

## Example Prompts

### Preranked GSEA on GO
> "I have a full DESeq2 result with the Wald statistic for every gene. Build a named decreasing ranked vector from the `stat` column, run GO biological-process GSEA with a fixed seed, and give me the top terms by adjusted p-value with their NES and leading-edge genes."

### ORA vs GSEA choice
> "I have a complete ranked DE result for all genes and no obvious significance cutoff. Decide whether ORA or GSEA is appropriate and run whichever fits on GO terms, explaining the choice."

### Ranking metric
> "Run GSEA using a signed p-value (`sign(log2FC) * -log10(p)`) as the ranking statistic instead of fold change, because my data came from edgeR which has no Wald-equivalent column."

### MSigDB and KEGG
> "Run GSEA against the MSigDB Hallmark collection on my ranked human gene list, then also run gseKEGG and note that KEGG queries the live database."

### Matrix and design available
> "I have a log-expression matrix and a case/control design. Test the KEGG pathways for a group difference in a way that is not fooled by co-regulated genes."

### Per-sample scores
> "Convert my expression matrix into a per-sample pathway-activity matrix with GSVA so I can cluster samples and correlate the scores with survival."

## Related Skills

- go-enrichment - Gene-list ORA alternative when no ranking exists
- kegg-pathways - KEGG pathway and module enrichment and ID conventions
- reactome-pathways - Reactome curated-pathway enrichment (local DB)
- wikipathways - WikiPathways community-pathway enrichment
- enrichment-visualization - gseaplot2, ridgeplot, and dotplot of GSEA results
- differential-expression/de-results - Source of the ranking statistic and the padj column conventions
- workflows/expression-to-pathways - End-to-end DE-to-enrichment pipeline
