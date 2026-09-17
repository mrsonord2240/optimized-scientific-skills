# Reactome Pathway Enrichment - Usage Guide

## Overview
ReactomePA tests a gene list (ORA via `enrichPathway`) or a ranked gene vector (GSEA via `gsePathway`) against Reactome, a curated, peer-reviewed knowledgebase whose atomic unit is the REACTION and whose pathways are nested containers of reactions. The result is therefore one signal projected onto a hierarchy, not a list of independent findings, and the honest read is the deepest curated human pathway nodes the list over-represents, deduplicated against their ancestors, against a measured background. ReactomePA reads the local `reactome.db`, so a run is reproducible offline given the Bioconductor release - unlike KEGG and WikiPathways, which query a live database.

Package install commands live in SKILL.md's Version Compatibility section; ID-prep, universe, organism-ceiling and hierarchy rules live in SKILL.md's Decision Tree and Per-Method Failure Modes sections.

## Example Prompts

### Over-representation
> "I have 240 significant genes as gene symbols from a DESeq2 contrast and the ~14,000 expressed genes as the background. Convert to Entrez, run Reactome over-representation, and give me the top pathways with fold enrichment, deduplicated so I am not double-counting parent and child pathways."

### GSEA
> "I have a full DESeq2 result with the test statistic for every gene. Build the ranked vector, fix the seed, and run Reactome GSEA, then show me the leading-edge genes for the top pathways."

### Hierarchy interpretation
> "My Reactome result has 'Cell Cycle Checkpoints', 'G2/M Checkpoints', and 'G1/S Transition' all near the top. Are those separate findings, and which one should I report?"

### Inspecting one pathway
> "Draw the reaction network for my top Reactome pathway colored by log2 fold change, and give me the link to open it in the Reactome Pathway Browser."

### Comparative / multi-omics
> "I have RNA-seq counts for treated vs control. Use ReactomeGSA to find which Reactome pathways differ between the groups."

## Reactome vs KEGG

| Feature | Reactome | KEGG |
|---------|----------|------|
| Atomic unit | Reaction (typed entities, PubMed-cited) | Pathway map |
| Curation | Expert-authored, externally peer-reviewed | KEGG-team curated |
| Structure | Deep event hierarchy (parent/child double-count) | Mostly flat map list |
| Granularity | Finest (reaction level); heavier multiple-testing | Coarser (pathway level) |
| Reproducibility | Local reactome.db, pinned to the Bioconductor release | Live REST API, date-dependent |
| Metabolic depth | Good | Deeper |
| License | CC0 (fully open) | Free academic; commercial license for KEGG REST |
| Organisms (in R) | 7 in ReactomePA (DB projects to ~14-20) | 8,000+ |

The results-column reference (`enrichResult`/`gseaResult`), failure modes (ENTREZ-only, universe, hierarchy, viewPathway, non-human orthology, R-vs-web skew, organism ceiling), and the ReactomeGSA comparative pointer all live in SKILL.md (Understanding Results, Per-Method Failure Modes, Common Errors). See enrichment-visualization for dotplot/emapplot/cnetplot/gseaplot2 recipes.

## Related Skills

- go-enrichment - The hypergeometric test and the background-universe problem
- gsea - The GSEA running-sum engine and ranking-metric choice
- kegg-pathways - KEGG pathway/module enrichment; deeper metabolic coverage
- wikipathways - WikiPathways community-pathway enrichment (also CC0)
- enrichment-visualization - Dot/bar/cnet/emap/tree/GSEA plots of enrichment results
- differential-expression/de-results - Source of the gene list and the ranking statistic
- workflows/expression-to-pathways - End-to-end DE-to-enrichment pipeline
