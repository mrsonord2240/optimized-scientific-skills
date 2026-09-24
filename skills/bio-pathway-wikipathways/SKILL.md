---
name: bio-pathway-wikipathways
category: Data Analysis
description: Tests a gene list (ORA, enrichWP) or a ranked gene vector (GSEA, gseWP) against the WikiPathways community-curated pathway collection with clusterProfiler and rWikiPathways. Covers why a WikiPathways result is a snapshot of a live, monthly-updated database (enrichWP/gseWP/gson_WP silently pull data.wikipathways.org/current/), why reproducibility requires pinning a dated GMT via downloadPathwayArchive(date=, format='gmt'), why the WP GMT is Entrez-keyed so symbols and Ensembl silently overlap nothing, why universe=NULL gives a biased all-WP-genes background, how to split the name%version%wpid%org term, and why WikiPathways (CC0, no peer review) complements KEGG/Reactome. Use when running open community-pathway enrichment, covering a non-model WP species, catching disease/drug pathways missing from KEGG/Reactome, or needing a reproducible dated analysis. The gene list comes from differential-expression/de-results; visualize with enrichment-visualization.
tool_type: r
primary_tool: rWikiPathways
license: MIT
author: GPTomics
---

## Version Compatibility

Reference examples tested with: clusterProfiler 4.18+, rWikiPathways 1.26+, org.Hs.eg.db 3.18+.

```r
if (!require('BiocManager', quietly = TRUE))
    install.packages('BiocManager')

BiocManager::install(c('clusterProfiler', 'rWikiPathways', 'enrichplot', 'org.Hs.eg.db'))
```

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

WikiPathways is a LIVE, monthly-updated database. `enrichWP`, `gseWP`, and `gson_WP` all download `data.wikipathways.org/current/gmt/` at run time, so the SAME code returns DIFFERENT pathways and p-values months apart with no error. `current/` is not a version. For a reproducible analysis pin a dated release with `downloadPathwayArchive(date='YYYYMMDD', organism=, format='gmt')` and report the date. All `enrichWP`/`gseWP`/`downloadPathwayArchive` calls require internet at run time.

# WikiPathways Enrichment

**"Which community-curated WikiPathways are enriched in my genes?"** -> Test WikiPathways gene sets against a gene list (ORA) or a ranked vector (GSEA), pinning a dated GMT for reproducibility - because the live monthly database changes under identical code, and the WP GMT is Entrez-keyed so any other ID type silently overlaps nothing.
- R (ORA): `enrichWP(entrez, organism='Homo sapiens', universe=all_entrez)`
- R (GSEA): `gseWP(named_decreasing_entrez_vector, organism='Homo sapiens')`
- R (reproducible): `downloadPathwayArchive(date='YYYYMMDD', organism=, format='gmt')` -> `read.gmt` -> split term -> `enricher`/`GSEA`

Scope: WikiPathways-specific enrichment - the data model, the `current/`-vs-dated GMT reproducibility pin, the Entrez-GMT requirement, the term-field split, PFOCR as a noisier complement, and the WP-vs-KEGG-vs-Reactome contrast. The ORA/GSEA method choice and hypergeometric/background theory -> the category README. The DE list and ranking statistic -> differential-expression/de-results. KEGG and Reactome -> kegg-pathways, reactome-pathways. Plot grammar -> enrichment-visualization.

## The Single Most Important Modern Insight -- A WikiPathways Result Is a Snapshot of a Live, Community-Edited Database Taken on the Run Date

WikiPathways is a wiki: anyone can create or edit a pathway, content is CC0, and there is NO formal journal-style peer review gating a pathway's publication (Pico 2008 *PLoS Biol* 6:e184; Martens 2021 *NAR* 49:D613). The collection is republished as a dated GMT archive every MONTH. Three properties every misuse forgets:

1. **`current/` is not a version.** `enrichWP`, `gseWP`, and `gson_WP` all silently download `data.wikipathways.org/current/gmt/` - the latest monthly release. Identical code two months apart returns different pathways and different p-values, with no error and no warning. Reproducibility is NOT a code freeze; it is a dated GMT: `downloadPathwayArchive(date='YYYYMM10', organism='Homo sapiens', format='gmt')` (releases land on the 10th of each month), read it, run `enricher`/`GSEA` on the pinned sets, and report the date in methods. **The live archive keeps only the last ~12 months** (data.wikipathways.org's own landing page: "This site hosts monthly data releases for the last 12 months") - never hardcode an old literal date, compute a recent one instead (e.g. `format(Sys.Date() - 60, '%Y%m10')`), and for anything older use the permanent, citable Zenodo GMT/GPML archive (https://zenodo.org/communities/wikipathways). `gson_WP()` freezes only within a session (it snapshots `current/`), not across time.

2. **The WP GMT speaks Entrez, and the wrong ID type fails silently.** The GMT is Entrez-keyed via BridgeDb. Passing SYMBOL or ENSEMBL yields near-zero overlap and an empty or misleading result with NO error - convert to Entrez upstream (`bitr`/OrgDb) before `enrichWP`. Likewise `universe=NULL` makes the background "all genes that happen to be in WP" - a small, biased set that inflates significance; pass the assayed/tested Entrez vector as `universe`.

3. **A community pathway is a hypothesis someone drew, not a reviewed fact.** The two things that make WP valuable (open CC0 license, anyone-can-edit curation that captures disease/drug pathways KEGG and Reactome lack, e.g. the COVID-19 Disease Map) are the same two things that make its quality heterogeneous. Many WP pathways are also imported from KEGG/Reactome, so "three databases agree" can be circular rather than independent. Treat each hit as a community claim - check `getPathwayInfo(WPID)` last-edit/curation before leaning on a single WP pathway for a key conclusion - and run WP as a COMPLEMENT to KEGG/Reactome, never a sole peer-reviewed source.

## Tool Taxonomy

| Source / function | Citation | Mechanism / role | When |
|-------------------|----------|------------------|------|
| WikiPathways ORA (`enrichWP`) | Pico 2008 *PLoS Biol* 6:e184; Martens 2021 *NAR* 49:D613; Wu 2021 *Innovation* 2:100141 | hypergeometric test vs the WP GMT (delegates to `enricher`); downloads `current/` | a thresholded gene LIST against community pathways |
| WikiPathways GSEA (`gseWP`) | Agrawal 2024 *NAR* 52:D679; Wu 2021 *Innovation* 2:100141 | running-sum FCS over a ranked vector (delegates to `GSEA`); downloads `current/` | all genes ranked, no arbitrary cutoff |
| `rWikiPathways` (query/download) | Slenter/Hanspers/Pico, Bioconductor | API client: `listOrganisms`, `listPathways`, `getPathwayInfo`, `getXrefList`, `findPathwaysByText`, `downloadPathwayArchive` | inspect pathways, fetch genes, pin a dated GMT |
| Dated GMT + `enricher`/`GSEA` | Wu 2021 *Innovation* 2:100141 | run enrichment on a pinned, parsed GMT, bypassing auto-download | the REPRODUCIBLE pattern; report the date |
| PFOCR (Pathway Figure OCR) | Hanspers 2020 *Genome Biol* 21:273; Shin 2023 *BMC Genomics* 24:713 | machine-OCR'd gene sets from published figures; larger + noisier, no edges | high-recall disease/process coverage as a complement; NOT what `enrichWP` queries |
| KEGG / Reactome (siblings) | -> kegg-pathways, reactome-pathways | metabolic/signaling maps (live) / curated reactions (local) | the primary databases WP complements |

## WikiPathways vs KEGG/Reactome

| Feature | WikiPathways | KEGG | Reactome |
|---------|--------------|------|----------|
| License | CC0 (fully open) | Restrictive (commercial bulk/API) | CC-BY / CC0 |
| Curation | Community wiki, no formal peer review | Largely automated KO reconstruction | Expert-curated and reviewed |
| Species | ~30+ | 4000+ (genome-derived) | ~15 (deep human) |
| Human pathways | ~1,100 (`listPathways('Homo sapiens')`, 2026-09) | ~370 (KEGG REST `list/pathway/hsa`, 2026-09) | see reactome-pathways |
| Focus | Disease/drug + general | Metabolic/signaling | Reaction-level mechanism |
| Reproducibility | Pin a dated monthly GMT (live `current/` otherwise) | Live REST API (date-dependent) | Local reactome.db (version-pinned) |

## Decision Tree by Scenario

| Scenario | Recommended | Why |
|----------|-------------|-----|
| Quick exploratory ORA, reproducibility not yet needed | `enrichWP(entrez, organism, universe=all_entrez)` | fastest path; log that it used the `current/` release |
| Publication / reproducible analysis | `downloadPathwayArchive(date='YYYYMMDD', organism, format='gmt')` -> read -> split -> `enricher`/`GSEA` | the dated GMT is the only cross-time pin; report the date |
| All genes carry a DE statistic, cutoff would be arbitrary | `gseWP` (see the category README for the ORA/GSEA choice) | ranked FCS uses the full list, no cutoff |
| Pre-selected list (module, screen hits, GWAS loci) | `enrichWP` ORA | no ranking available |
| Disease / drug pathways missing from KEGG/Reactome | WP as a complement, run alongside KEGG/Reactome | community content is genuinely additive where it exists |
| Maximum gene/process coverage, noise tolerable | PFOCR (separate resource), not `enrichWP` | figure-OCR sets are higher-recall, lower-precision |
| Non-model but WP-supported species (zebrafish, fly, worm, Arabidopsis) | `enrichWP(entrez, '<scientific name>')`, verify via `rWikiPathways::listOrganisms()` | WP covers ~30+ species |
| Compare up- vs down-regulated | `compareCluster(geneClusters=list(up=..,down=..), fun='enrichWP', organism=)` | one model, faceted dotplot |
| Genes are SYMBOL/ENSEMBL | convert to Entrez first (`bitr`) | the WP GMT is Entrez-keyed; other types overlap nothing |

## Agent Workflow

1. Load DE results and extract the significant gene list (ORA) or build the named decreasing ranking vector (GSEA).
2. Convert gene IDs to Entrez with `bitr` and set the background universe to the tested genes.
3. For a reproducible run, pin a dated GMT with `downloadPathwayArchive(date=, format='gmt')`, split the `name%version%wpid%org` term field, and run `enricher`/`GSEA`; otherwise run `enrichWP`/`gseWP` and log that it used the `current/` release.
4. Make the result readable with `setReadable()` and report `p.adjust`/`qvalue` (not raw p).
5. Hand the result object to enrichment-visualization for plots and corroborate hits against KEGG/Reactome.

## Over-Representation Analysis (enrichWP)

**Goal:** Find WikiPathways over-represented in a thresholded gene list, against a defensible background.

**Approach:** Convert significant genes to Entrez, pass the tested-gene set as `universe`, run `enrichWP`, then make the result readable. `enrichWP` downloads the `current/` GMT - acceptable for exploration, but pin a date for anything reportable.

```r
library(clusterProfiler)
library(org.Hs.eg.db)

# enrichWP downloads the current/ WP GMT over the network; symbols/Ensembl must be Entrez first
sig <- bitr(sig_symbols, fromType='SYMBOL', toType='ENTREZID', OrgDb=org.Hs.eg.db)$ENTREZID
all_entrez <- bitr(tested_symbols, fromType='SYMBOL', toType='ENTREZID', OrgDb=org.Hs.eg.db)$ENTREZID

wp <- enrichWP(gene=sig, organism='Homo sapiens', universe=all_entrez,
               pvalueCutoff=0.05, pAdjustMethod='BH', minGSSize=10, maxGSSize=500, qvalueCutoff=0.2)
wp <- setReadable(wp, OrgDb=org.Hs.eg.db, keyType='ENTREZID')   # geneID column -> symbols
as.data.frame(wp)   # ID=WPID, Description, GeneRatio, BgRatio, p.adjust, qvalue, Count
```

## GSEA (gseWP)

**Goal:** Find WikiPathways whose genes shift coordinately across the full ranking, with no cutoff.

**Approach:** Build a NAMED Entrez vector sorted DECREASING by the ranking metric, fix the permutation seed, then run `gseWP`. There is no `universe` argument - FCS uses the whole ranked list.

```r
gl <- sort(setNames(de$log2FoldChange, de$entrez), decreasing=TRUE)   # named, decreasing, Entrez names
set.seed(123)                                                          # fix permutation reproducibility
wp_gsea <- gseWP(geneList=gl, organism='Homo sapiens',
                 pvalueCutoff=0.05, pAdjustMethod='BH', minGSSize=10, maxGSSize=500)
as.data.frame(wp_gsea)   # NES, p.adjust, core_enrichment (the leading edge)
```

## Reproducible Analysis with a Dated GMT (the correct pattern)

**Goal:** Make a WP analysis reproducible across re-runs by pinning a dated release instead of pulling `current/`.

**Approach:** Download a dated GMT (pass `format='gmt'` - the default is `gpml`), split the compound `name%version%wpid%org` term field into TERM2GENE/TERM2NAME, run `enricher`/`GSEA` on the pinned sets, and report the date in methods. The live archive retains only the last ~12 months of monthly releases (10th of each month) - compute a recent date rather than hardcoding one that will 404 as time passes, and loop over successive months (`Sys.Date() - 60, -90, -120, ...`) with `tryCatch` so a single missing release does not stop the run - `scripts/wikipathways_pinned_enrich.R` does this and prints the date that succeeded. If every in-window date fails, or for a fixed historical date beyond the window, use the Zenodo GMT/GPML archive instead (https://zenodo.org/communities/wikipathways).

```bash
# sig_entrez.txt / universe_entrez.txt: one Entrez ID per line (tested genes as universe); prints the release date used
Rscript scripts/wikipathways_pinned_enrich.R sig_entrez.txt universe_entrez.txt 'Homo sapiens' wp_pinned.csv
```

For GSEA, build the same `t2g` (`wpid`, `gene`) / `t2n` (`wpid`, `name`) tables and call `GSEA(geneList, TERM2GENE=t2g, TERM2NAME=t2n)`; `examples/wikipathways_explore.R` shows the pattern inline.

`gson_WP(organism)` returns a GSON snapshot object, but it still pulls `current/` - it freezes a session, NOT a chosen historical date. Only the dated `downloadPathwayArchive` GMT survives a re-run months later.

## Query the Database Directly (rWikiPathways)

```r
library(rWikiPathways)

organisms <- rWikiPathways::listOrganisms()  # supported species (full scientific names; ~30+)
organisms
listPathways('Homo sapiens')             # all WPIDs + names for a species
getPathwayInfo('WP554')                  # metadata incl. last-edit; check before trusting a single hit
getXrefList('WP554', 'L')                # genes by BridgeDb system code: 'L'=Entrez, 'H'=HGNC, 'En'=Ensembl
findPathwaysByText('cancer')             # text search (searchPathways() is NOT a current function)
```

## Other Organisms

```r
wp_mouse <- enrichWP(gene=mouse_entrez, organism='Mus musculus')
wp_zfish <- enrichWP(gene=zfish_entrez, organism='Danio rerio')
# verify the exact organism string before running:
rWikiPathways::listOrganisms()
```

## Understanding Results

| Column | Description |
|--------|-------------|
| ID | WikiPathways stable ID (WP####) |
| Description | Pathway name |
| GeneRatio | Query genes in the pathway / query genes mapped to any pathway |
| BgRatio | Pathway genes in the universe / universe genes mapped |
| pvalue | Raw p-value |
| p.adjust | BH-adjusted p-value (report this, not raw p) |
| qvalue | q-value |
| geneID | Genes in the pathway (symbols after setReadable) |
| Count | Number of query genes in the pathway |

For GSEA results read `NES` (sign = direction along the ranking) and `core_enrichment` (the leading-edge genes).

## Per-Method Failure Modes

### Unpinned current/ release
**Trigger:** running `enrichWP`/`gseWP`/`gson_WP` without `downloadPathwayArchive(date=)`. **Mechanism:** all three download `data.wikipathways.org/current/`, the latest monthly release. **Symptom:** the same script returns different pathways/p-values months apart, with no error. **Fix:** pin a dated GMT, run `enricher`/`GSEA` on it, and report the date.

### Symbols or Ensembl into an Entrez GMT
**Trigger:** passing SYMBOL/ENSEMBL IDs to `enrichWP`/`gseWP`. **Mechanism:** the WP GMT is Entrez-keyed via BridgeDb, so non-Entrez IDs overlap nothing. **Symptom:** an empty or near-empty result, NO error. **Fix:** `bitr` to ENTREZID first; confirm the conversion rate before trusting the result.

### Default universe inflates significance
**Trigger:** `universe=NULL` (the default). **Mechanism:** `enricher` then uses "all genes in the WP GMT" as background - a small, biased set, not the assayed genes. **Symptom:** implausibly strong p-values for tissue-specific or off-target pathways. **Fix:** pass the tested-gene Entrez vector as `universe`.

### gson_WP mistaken for a reproducibility pin
**Trigger:** treating `gson_WP()` as "the snapshot" for a reproducible analysis. **Mechanism:** it snapshots `current/` into an object - it freezes a session, not a historical date. **Symptom:** a re-run months later gives a different snapshot. **Fix:** use the dated `downloadPathwayArchive` GMT for cross-time reproducibility.

### Unsplit GMT term field
**Trigger:** `read.gmt` on a WP GMT without splitting the term. **Mechanism:** the set-name field is a compound `name%version%wpid%org` joined by `%`. **Symptom:** WPIDs and clean names are buried in one column; TERM2GENE/TERM2NAME are wrong. **Fix:** `separate(., term, c('name','version','wpid','org'), sep='%')` (or use `read.gmt.wp`).

### searchPathways() is gone
**Trigger:** calling `searchPathways('cancer', 'Homo sapiens')`. **Mechanism:** it is not a current rWikiPathways function. **Symptom:** an error. **Fix:** `findPathwaysByText()` / `findPathwayIdsByText()`.

### format defaults to gpml
**Trigger:** `downloadPathwayArchive(date=, organism=)` without `format='gmt'`. **Mechanism:** `format` defaults to `gpml`, which `read.gmt` cannot read. **Symptom:** a GPML file or a parse error. **Fix:** pass `format='gmt'`.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| `pvalueCutoff=0.05` | `enricher`/`GSEA` default | filters on p.adjust (BH) by default; standard FDR gate |
| `qvalueCutoff=0.2` | clusterProfiler `enricher` default | secondary q-value gate on ORA |
| `pAdjustMethod='BH'` | clusterProfiler default | Benjamini-Hochberg FDR; not Bonferroni (too conservative for gene-set screens) |
| `minGSSize=10` | `enricher`/`GSEA` default | drop tiny WP pathways that overfit; many WP specialist sets fall below this and are never tested |
| `maxGSSize=500` | `enricher`/`GSEA` default | drop overly broad sets that always "enrich" |
| `set.seed(123)` for `gseWP` | reproducibility convention | permutation p-values drift across runs without a fixed seed (any fixed seed works) |
| Pin `date='YYYYMMDD'` | Martens 2021 *NAR* 49:D613 | WP republishes monthly; `current/` is not a version, so report the dated release |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| `enrichWP` returns 0 terms | passed SYMBOL/ENSEMBL not Entrez | `bitr` to ENTREZID first |
| Implausibly significant pathways | `universe=NULL` (all-WP-genes background) | pass the tested-gene Entrez vector as `universe` |
| Different results each run | unpinned `current/` release | `downloadPathwayArchive(date=, format='gmt')`; report the date |
| `searchPathways` error | function removed | use `findPathwaysByText()` |
| `read.gmt` term column is a `%`-compound | term field not split | `separate(., term, c('name','version','wpid','org'), sep='%')` |
| `downloadPathwayArchive` opens a browser / downloads nothing | `organism=NULL` | name the organism to actually download a file |
| `downloadPathwayArchive` warns `cannot open URL ... 404` then errors | that month's release is outside the ~12-month window or was not published | step back one release (`Sys.Date() - 90, -120, ...`, 10th of the month); for an older fixed date use the Zenodo archive |
| GPML where a GMT was expected | `format` defaulted to `gpml` | pass `format='gmt'` |
| `gseWP` error about vector names | geneList not named or not sorted decreasing | build a named Entrez vector, `sort(decreasing=TRUE)` |
| `enrichWP`/`gseWP` returns NULL with no terms | wrong/non-canonical organism string (e.g. a common name like `'zebrafish'`) | verify with `rWikiPathways::listOrganisms()` first; the string must match exactly |

## References

- Pico AR, Kelder T, van Iersel MP, Hanspers K, Conklin BR, Evelo C. 2008. WikiPathways: pathway editing for the people. *PLoS Biol* 6(7):e184.
- Martens M, Ammar A, Riutta A, et al. 2021. WikiPathways: connecting communities. *Nucleic Acids Res* 49(D1):D613-D621.
- Agrawal A, Balci H, Hanspers K, et al. 2024. WikiPathways 2024: next generation pathway database. *Nucleic Acids Res* 52(D1):D679-D689.
- Hanspers K, Riutta A, Summer-Kutmon M, Pico AR. 2020. Pathway information extracted from 25 years of pathway figures. *Genome Biol* 21:273.
- Shin MG, Pico AR. 2023. Using published pathway figures in enrichment analysis and machine learning. *BMC Genomics* 24:713.
- Wu T, Hu E, Xu S, et al. 2021. clusterProfiler 4.0: A universal enrichment tool for interpreting omics data. *The Innovation* 2(3):100141.

## Related Skills

- go-enrichment - GO over-representation alternative
- gsea - Ranked-list GSEA mechanics and the ranking metric
- kegg-pathways - KEGG pathway/module enrichment (the primary DB WP complements)
- reactome-pathways - Reactome curated-pathway enrichment (the primary DB WP complements)
- enrichment-visualization - Dot/bar/cnet/emap/GSEA plots of the enrichment result
- differential-expression/de-results - Source of the gene list and the ranking statistic
- workflows/expression-to-pathways - End-to-end DE-to-enrichment pipeline
