# KEGG Pathway and Topology Enrichment - Usage Guide

## Overview
KEGG enrichment tests genes against KEGG's curated pathway and module gene sets across all three generations of pathway analysis: over-representation (enrichKEGG, enrichMKEGG), ranked GSEA (gseKEGG), and signed-topology perturbation (SPIA/graphite). KEGG is the database that owns the third generation because it ships signed directed signaling topology (KGML), letting SPIA propagate fold-changes through the wiring rather than treating a pathway as an unordered gene set. The caveat that breaks reproducibility: a KEGG result is a timestamped query against a live, partially-paywalled REST API, so it is irreproducible unless the release is pinned.

Package install commands live in SKILL.md's Version Compatibility section; ID-prep, universe and pinning rules live in SKILL.md's Decision Tree, Agent Workflow, and Prepare the Gene IDs sections.

## Example Prompts

### KEGG over-representation
> "I have 240 significant genes from a DESeq2 contrast as SYMBOLs and the ~13,000 expressed genes as background. Convert both to Entrez, run KEGG pathway ORA for human with the measured universe, and give me the top pathways by adjusted p-value with fold enrichment."

### Signed-topology perturbation
> "I have a human DE list with log2 fold-changes and a universe. Run SPIA so direction and network position are used, and tell me which signaling pathways are activated vs inhibited - and explain why this is not appropriate for metabolic pathways."

### Prokaryotic / non-model
> "This is a Pseudomonas aeruginosa RNA-seq DE list with PA-locus-tag gene IDs. Run KEGG enrichment with the right organism code and keyType, without forcing an OrgDb or bitr."

### Reproducibility
> "Pin the current human KEGG release as a snapshot, record the date, and run my enrichment against the snapshot so a rerun next year gives the same pathways."

### Multi-condition and modules
> "Compare KEGG enrichment between my up- and down-regulated gene sets in one faceted dotplot, and also run KEGG module enrichment to localize which sub-process is hit."

The step-by-step agent workflow (including when to warn about a missing background gene set), organism-code table, and results-column reference all live in SKILL.md (Agent Workflow, Common Organism Codes, Understanding Results).

## Related Skills

- go-enrichment - Hypergeometric ORA and the background-universe problem
- gsea - GSEA running-sum engine and ranking-metric choice
- reactome-pathways - Reactome curated-pathway enrichment (reproducible local DB)
- wikipathways - WikiPathways community-pathway enrichment
- enrichment-visualization - Dot/bar/cnet/emap plots of enrichment results
- differential-expression/de-results - Source of the gene list and the fold-changes
- workflows/expression-to-pathways - End-to-end DE-to-enrichment pipeline
