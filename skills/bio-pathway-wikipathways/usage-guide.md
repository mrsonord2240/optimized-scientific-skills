# WikiPathways Enrichment - Usage Guide

## Overview
WikiPathways is an open, community-curated pathway database (CC0 license, ~30+ species, no formal peer review) built on a wiki model. This skill runs over-representation analysis (`enrichWP`) and GSEA (`gseWP`) against it with clusterProfiler, and uses rWikiPathways to query the database and pin a dated GMT. The point that governs reproducibility: a WikiPathways result is a snapshot of a live, monthly-updated database - the same code returns different results months apart unless a dated release is pinned. WikiPathways is a complement to KEGG/Reactome, not a sole source.

Package install commands live in SKILL.md's Version Compatibility section; ID-prep, universe, and dated-GMT pinning rules live in SKILL.md's Over-Representation Analysis, GSEA, and Reproducible Analysis sections.

## Example Prompts

### Basic enrichment
> "I have ~200 significant human genes as symbols and the ~13,000 tested genes as the background. Convert them to Entrez, run WikiPathways over-representation with the tested set as the universe, and give me the top 15 pathways by adjusted p-value with the gene symbols readable."

### Reproducible analysis
> "Run WikiPathways enrichment but make it reproducible: pin a dated GMT release, split the term field into the WPID and name, run enrichment on the pinned sets, and tell me which release date to report in the methods."

### ORA vs GSEA
> "I have a full ranked DESeq2 result for every gene with no clear cutoff. Should I run WikiPathways ORA or GSEA, and run whichever is appropriate."

### Organism-specific
> "Run WikiPathways enrichment for zebrafish Entrez genes, and first confirm the exact organism string WikiPathways expects."

### Combining databases
> "Run enrichment against WikiPathways, KEGG, and Reactome and tell me which pathways are unique to WikiPathways versus shared."

### Exploring pathways
> "Search WikiPathways for cancer-related pathways and show me the last-edited date for the top hit before I trust it."

The step-by-step agent workflow, the WikiPathways-vs-KEGG/Reactome comparison, and the results-column reference all live in SKILL.md (Agent Workflow, WikiPathways vs KEGG/Reactome, Understanding Results).

## Related Skills
- go-enrichment - Gene Ontology enrichment
- kegg-pathways - KEGG pathway enrichment
- reactome-pathways - Reactome pathway enrichment
- gsea - Gene Set Enrichment Analysis mechanics
- enrichment-visualization - Visualization of enrichment results
- differential-expression/de-results - Source of the gene list and ranking statistic
- workflows/expression-to-pathways - End-to-end DE-to-enrichment pipeline
