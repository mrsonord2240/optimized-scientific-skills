---
name: ca-bio-pathway-go-enrichment
category: Data Analysis
description: Use when finding which Gene Ontology terms are over-represented in a pre-selected gene list such as DE hits, a co-expression module or screen hits.
tool_type: r
primary_tool: clusterProfiler
license: MIT
author: GPTomics
---

# GO Over-Representation Analysis

Pick the one row that matches your data. Read that file and only that file. Run its command before writing any code of your own.
Paths are relative to this Skill's directory.

| Your data and task | Read |
|--------------------|------|
| A gene list plus the list of genes that could have been selected | `routes/ora.md` |
| A DE results table (padj, log2FoldChange, pvalue) | `routes/de-table.md` |
| RNA-seq DE hits and you want gene-length bias corrected (human, Ensembl IDs) | `routes/length-bias.md` |
| Your own gene sets (MSigDB subset, in-house) instead of GO | `routes/custom-sets.md` |
| A descriptive GO-slim breakdown of a list, no test | `routes/groupgo.md` |
| A ranked list of all genes with no cutoff, KEGG, Reactome, WikiPathways, or plots | `routes/other-skills.md` |
| Empty table, wrong ontology, "no gene can be mapped", redundant terms | `routes/errors.md` |

Rules on every route:

- The universe is the null. Always supply the genes that could have entered the list (tested or detected genes), mapped the same way as the foreground. Never run without one, and never use the whole genome unless every gene could have been detected.
- Set the ontology explicitly (`enrichGO` defaults to MF, not BP).
- `pvalueCutoff` filters `p.adjust`, not the raw p. Judge on `p.adjust` and read `FoldEnrichment` beside it.
- Run up- and down-regulated genes separately when direction matters.
- Collapse redundant ancestor terms (`simplify`) within one ontology; never count a lineage as independent hits.
- ORA gives hypotheses. It cannot validate the DE list that produced it.
- Record the OrgDb and GO.db versions with the results (the scripts write them to `results.csv.log`).

Tested with clusterProfiler 4.18.4, org.Hs.eg.db 3.22 (GOseq route: Bioconductor 3.20 pins goseq 1.58.0).
