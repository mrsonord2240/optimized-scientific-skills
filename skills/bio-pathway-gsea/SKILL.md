---
name: bio-pathway-gsea
description: Tests a ranked gene vector for coordinated expression shifts in GO, KEGG, Reactome, or MSigDB gene sets with clusterProfiler's gseGO, gseKEGG, gsePathway, and GSEA (fgseaMultilevel engine), and scores per-sample pathway activity with ssGSEA and GSVA. Covers why a GSEA result is a deterministic function of three implicit choices (the ranking STATISTIC, the weight exponent p, and which LABELS are permuted), why the input must be a NAMED vector sorted DECREASING by a signed variance-calibrated metric (DESeq2 stat, limma t) not a raw p-value that erases direction, why preranked gene-permutation is anti-conservative for correlated sets (CAMERA with inter.gene.cor=NA is the fix), why nPerm silently downgrades the engine instead of erroring (eps governs tiny p), and why set.seed is required. Use when every gene carries a DE statistic, when a hard cutoff is arbitrary, or when ORA finds nothing. For gene-list ORA see go-enrichment; the ranking statistic comes from differential-expression/de-results.
tool_type: r
primary_tool: clusterProfiler
license: MIT
author: GPTomics
---

## Version Compatibility

Verified 2026-09-23 in an isolated Linux environment: R 4.5.3 / Bioconductor 3.22, clusterProfiler 4.18.4, fgsea 1.36.2, org.Hs.eg.db 3.20.0, msigdbr 26.1.1, limma 3.66.0, GSVA 2.4.9, ReactomePA 1.54.0, and reactome.db 1.95.0. The shipped GO and MSigDB examples and representative GSVA, ID, KEGG, and Reactome calls must exit 0; materialized output with a nonzero process exit is not a successful run.

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

For a reproducible fresh runtime, create an isolated environment (do not repair a shared R library in place): `micromamba create -n gsea-r45 -c conda-forge -c bioconda r-base=4.5 bioconductor-clusterprofiler bioconductor-org.hs.eg.db bioconductor-fgsea bioconductor-limma bioconductor-reactomepa bioconductor-gsva r-msigdbr`. Alternatively install with `BiocManager::install(c('clusterProfiler','org.Hs.eg.db','fgsea','limma','ReactomePA','GSVA'))` and `install.packages('msigdbr')` in an isolated R library. Native Windows R 4.4.3 / Bioconductor 3.20 in the prior runtime is not supported for this Skill: package loads could complete the calculation but exit with code 2816.

gseKEGG queries the live KEGG REST API, so the same code returns different results as KEGG updates; pin the run date. gseGO/gsePathway and MSigDB GSEA use local annotation (`org.*.eg.db`, `reactome.db`, msigdbr) and are reproducible given the package version. The single source of truth for versions is this block, not headings.

# Gene Set Enrichment Analysis (GSEA)

**"Which pathways shift coordinately across my full ranked gene list, with no cutoff?"** -> Walk a weighted running-sum down the genome-wide ranking and test whether each gene set piles up at one END - because that score reports the structure of YOUR ranking, so the ranking metric and the permutation type, not the gene sets, decide the result.
- R: `gseGO(geneList, OrgDb, ont)`, `gseKEGG(geneList, organism)`, `GSEA(geneList, TERM2GENE)`

Scope: threshold-free Functional Class Scoring (FCS) of a RANKED vector - the running-sum ES, the ranking-metric choice, the permutation null, NES/FDR, the leading edge, the matrix-based CAMERA/fry tests, and per-sample ssGSEA/GSVA scores. A pre-selected unranked gene LIST -> go-enrichment (ORA). The ranking statistic source -> differential-expression/de-results. KEGG/Reactome/WikiPathways database semantics -> kegg-pathways, reactome-pathways, wikipathways. Plots -> enrichment-visualization.

## The Single Most Important Modern Insight -- A GSEA Result Is a Deterministic Function of Three Usually-Implicit Choices: the Ranking Statistic, the Weight Exponent p, and the Permuted Labels

GSEA is not a discovery about biology - it is the running-sum's report on the ranking it was handed. The weighted enrichment score (Subramanian 2005; the default `exponent=1` weights each hit by the gene's statistic magnitude) asks exactly one question: do the members of a set pile up at one END of YOUR ranking. So the result is fixed by three decisions tutorials usually leave silent, and Wijesooriya 2022 found most published GSEA papers report none of them.

1. **The ranking statistic IS the experiment.** Rank by a signed, variance-calibrated metric (DESeq2 `stat`, limma moderated `t`) - sign gives direction, variance-calibration sinks noisy low-information genes to the middle. Ranking by a RAW p-value erases the sign, so up- and down-regulated genes collapse together and NES becomes uninterpretable. Ranking by bare log2FC lets a handful of low-count genes with huge unstable fold changes hijack the leading edge. A bad ranking is faithfully reported as a ranking artifact.
2. **The permutation type sets validity.** Phenotype (sample-label) permutation preserves the gene-gene correlation that co-regulated pathways are made of - it is the gold standard for type-I error control because it preserves that correlation, but needs the expression matrix and adequate n per group (~>=7); below that n its validity degrades. clusterProfiler/fgsea preranked are FORCED into GENE permutation, which treats genes as independent, destroys that correlation, and is ANTI-CONSERVATIVE: it returns "significant" pathways that are nothing but co-expression. Accept it and report it when only a ranked vector exists; when the expression matrix and design are in hand, run CAMERA with `inter.gene.cor = NA` as the DEFAULT test instead (see "Matrix + Design"). Measured 2026-09-21 on 24 samples with 20 co-regulated KEGG sets carrying no group effect: preranked GSEA on limma `t` called **9/20** at FDR<0.05, CAMERA `inter.gene.cor=NA` **0/20**, CAMERA default preset 7/20, GSVA + limma on the scores 1/20 (planted sets recovered by all). Full competitive/self-contained theory is in the category README.
3. **The honesty is bounded by reporting.** Log the ranking metric, the exponent p, the permutation type, the gene-set collection and its version/date, the size filters, and the multiple-testing method. Without those, the result is unfalsifiable.

## Tool Taxonomy

| Method | Citation | Mechanism / role | When |
|--------|----------|------------------|------|
| Preranked GSEA (gseGO/gseKEGG/GSEA) | Subramanian 2005 *PNAS* 102:15545; Mootha 2003 *Nat Genet* 34:267 | weighted running-sum ES over a ranked vector; gene-permutation null | the common case: a ranked statistic for all genes, no matrix |
| fgsea engine | Korotkevich 2021 *bioRxiv* 060012 (preprint) | fgseaMultilevel; resolves tiny p accurately down to `eps` | the engine under `by='fgsea'` (default); what gives sub-1/nperm p-values |
| Phenotype-permutation GSEA | Subramanian 2005 *PNAS* 102:15545 | shuffles sample labels; preserves gene-gene correlation | matrix + phenotype + ~>=7/group; the gold-standard competitive test. Concept only: it needs the Broad GSEA desktop/CLI, which this Skill does not run - use CAMERA below |
| CAMERA | Wu & Smyth 2012 *NAR* 40:e133 | competitive, VIF-corrects inter-gene correlation analytically - only with `inter.gene.cor = NA`; the default preset (0.01) gives no protection | matrix + design; want a correlation-honest competitive test |
| ROAST / fry | Wu 2010 *Bioinformatics* 26:2176 | self-contained rotation test; valid at any n | matrix + design, tiny n, "is the set DE at all" |
| ssGSEA | Barbie 2009 *Nature* 462:108 | per-sample rank-based enrichment score | a per-sample pathway-activity matrix (no contrast test) |
| GSVA | Hanzelmann 2013 *BMC Bioinformatics* 14:7 | unsupervised per-sample, per-set kernel/CDF score | per-sample features for clustering/survival/ML |

## Decision Tree by Scenario

| Scenario | Recommended | Why |
|----------|-------------|-----|
| Matrix + design available, want a set-level contrast (DEFAULT) | `limma::camera(..., inter.gene.cor = NA)` - "Matrix + Design" | correlation-honest: 0/20 co-regulated nulls called vs 9/20 for preranked GSEA on the same data |
| Only a ranked statistic for ALL genes (no matrix), cutoff would be arbitrary | preranked `gseGO`/`gseKEGG`/`GSEA` | threshold-free, but gene-permutation is anti-conservative for co-regulated sets: report it as screening and say so |
| Pre-selected unranked list (module, GWAS hits, screen) | ORA -> go-enrichment | no genome-wide ranking exists; ORA finds strong individual changes, GSEA coordinated subtle shifts |
| Function annotation, broad GO coverage | `gseGO` | GO is the broadest local resource |
| Metabolic / signaling pathways | `gseKEGG` -> kegg-pathways | KEGG maps (live DB) |
| Reaction-level, reproducible offline | `gsePathway` -> reactome-pathways | local reactome.db |
| Curated MSigDB hallmark / C2 / C5 | `GSEA(TERM2GENE)` + msigdbr | generic-input GSEA on any collection |
| Matrix + design, tiny n / covariates, "is set DE at all" | `limma::fry`/`roast` - "Matrix + Design" | self-contained rotation, valid at any n |
| Per-sample pathway-activity matrix for clustering/ML | ssGSEA / GSVA | scores each sample; a contrast on the scores is limma on the score matrix |
| The DE statistic / ranking itself | -> differential-expression/de-results | upstream, not enrichment |

## Build the Ranked Vector (the Caller-Owned Step)

**Goal:** Turn a DE table into the named numeric vector sorted strictly decreasing that every preranked function requires, ranked by a signed variance-calibrated metric.

**Approach:** Pick the ranking metric to match the DE tool, name the vector by gene ID, drop NAs, deduplicate to one statistic per gene, and sort decreasing. The DE mechanics and the `$padj`/`$adj.P.Val` column conventions live in differential-expression/de-results.

| DE source | Ranking metric | Column | Why |
|-----------|----------------|--------|-----|
| DESeq2 | Wald statistic | `stat` | signed + variance-calibrated; best single choice for RNA-seq |
| limma / voom | moderated t-statistic | `t` | empirical-Bayes shrinkage borrows variance; signed |
| edgeR (QL) | `sign(logFC) * -log10(PValue)` | derived | no single signed statistic column |
| any tool, last resort | log2 fold change | `log2FoldChange` | magnitude only; noisy for low-count genes |

```r
library(clusterProfiler)
library(org.Hs.eg.db)

de <- read.csv('de_results.csv')          # DE list source: differential-expression/de-results
gene_list <- de$stat                       # DESeq2 Wald stat: signed + variance-calibrated
names(gene_list) <- de$entrez_id
gene_list <- gene_list[!is.na(gene_list)]
gene_list <- gene_list[!duplicated(names(gene_list))]   # one statistic per gene; duplicates double-count hits
gene_list <- sort(gene_list, decreasing = TRUE)         # REQUIRED: unsorted input silently mis-ranks
# limma instead: tt <- topTable(eBayes(lmFit(expr, design)), coef = 2, number = Inf, sort.by = 'none'); gene_list <- setNames(tt$t, rownames(tt))
```

Ranking by `sign(log2FC) * -log10(pmax(pvalue, 1e-30))` (for edgeR, or when a Wald stat is unavailable) preserves direction and clamps `p==0` from going to `Inf` - **the clamp constant is a ranking weight, not just an `Inf` guard** (numbers in Quantitative Thresholds; `1e-300` lost a planted down-regulated set that `1e-30` recovered). Never rank by raw p-value alone or by `lfcShrink(type='normal')` (failure modes below). apeglm/ashr-shrunk results DROP the `stat` column - pull `stat` from the unshrunk `results(dds)` if ranking by it.

## Run Preranked GSEA on GO

**Goal:** Find GO terms whose members shift coordinately up or down across the full ranking, with no significance cutoff.

**Approach:** Set the permutation seed for reproducibility, set `eps=0` for exact tiny p-values, run gseGO, then map the leading-edge IDs back to symbols and read `core_enrichment` as the interpretable core.

```r
set.seed(123)                              # fixes the multilevel Monte Carlo; any fixed seed
gse_go <- gseGO(geneList = gene_list, OrgDb = org.Hs.eg.db, keyType = 'ENTREZID',
                ont = 'BP', exponent = 1, minGSSize = 10, maxGSSize = 500,
                eps = 0, pvalueCutoff = 0.05, pAdjustMethod = 'BH',
                seed = TRUE, by = 'fgsea', verbose = FALSE)
gse_go <- setReadable(gse_go, OrgDb = org.Hs.eg.db, keyType = 'ENTREZID')
```

`gseaResult` columns: `ID`, `Description`, `setSize`, `enrichmentScore` (raw ES), `NES`, `pvalue`, `p.adjust` (BH), `qvalue`, `rank` (ES-peak position), `leading_edge`, `core_enrichment` (`/`-separated leading-edge IDs). Report `p.adjust`/`qvalue`, never raw `pvalue` and never an invented `$FDR`. NES sign = direction: positive = top of the ranking (up in the contrast), negative = bottom.

## Run GSEA on KEGG, Reactome, or MSigDB

**Goal:** Apply the same preranked engine to a pathway database with the gene-ID type that database expects.

**Approach:** gseKEGG/gsePathway need ENTREZ-style IDs; choose the collection, keep the seed and `eps=0`, and note KEGG queries the live REST API (date-dependent) while Reactome and MSigDB are local.

```r
set.seed(123)
gse_kegg <- gseKEGG(geneList = gene_list, organism = 'hsa', keyType = 'ncbi-geneid',   # Entrez-named vector; 'kegg' keyType is the prokaryote locus-tag path
                    minGSSize = 10, maxGSSize = 500, eps = 0,
                    pvalueCutoff = 0.05, seed = TRUE, verbose = FALSE)   # live KEGG API; pin the date

library(msigdbr)
h <- msigdbr(species = 'Homo sapiens', collection = 'H')    # 26.x: collection= (was category=); gs_collection (was gs_cat)
t2g <- h[, c('gs_name', 'ncbi_gene')]                       # 26.x Entrez column is ncbi_gene; older releases used entrez_gene
gse_h <- GSEA(geneList = gene_list, TERM2GENE = t2g, exponent = 1,
              minGSSize = 10, maxGSSize = 500, eps = 0,
              pvalueCutoff = 0.05, seed = TRUE, verbose = FALSE)
```

**The TERM2GENE column must match the ID type the ranked vector is named by:** `ncbi_gene` for an Entrez-named vector (above), `gene_symbol` when the vector is named by symbols (e.g. limma `t` from a symbol-rowname matrix): `t2g <- unique(h[, c('gs_name', 'gene_symbol')])`. A mismatch fails with `--> No gene can be mapped` then `'organism' is not a slot in class "NULL"` (checked on clusterProfiler 4.14.6). For KEGG use `msigdbr(species = 'Homo sapiens', collection = 'C2', subcollection = 'CP:KEGG_LEGACY')` (msigdbr 26.x; `subcollection = 'KEGG'` is unknown there).

If the installed msigdbr still uses `category=`/`entrez_gene`, the old form works but warns - check `?msigdbr` and `names(h)`. ReactomePA's `gsePathway(geneList, organism='human')` reads the local reactome.db and also needs ENTREZ IDs.

## Matrix + Design: CAMERA, fry, and Contrasts on Scores

**Goal:** Test each gene set for a case-vs-control shift when the expression matrix and design are available, without the gene-independence assumption that makes preranked GSEA anti-conservative.

**Approach:** Map the gene sets to row indices, then run `camera` with the inter-gene correlation estimated (`inter.gene.cor = NA`). `fry` is the self-contained alternative for tiny n or "is the set DE at all". Checked on limma 3.62.2 with 24 samples x 8,288 genes and 186 KEGG sets: CAMERA `NA` recovered the planted up set (FDR 3e-17) and down set (FDR 2e-17) and called 0 of 20 co-regulated null sets; `fry` called 1 of 20.

```r
library(limma)
library(msigdbr)
# expr_matrix: log-scale genes x samples, rownames = symbols; group: factor(Control, Case)
design <- model.matrix(~ group)
kegg <- msigdbr(species = 'Homo sapiens', collection = 'C2', subcollection = 'CP:KEGG_LEGACY')
gene_sets <- lapply(split(kegg$gene_symbol, kegg$gs_name), unique)   # named list of symbol vectors
idx <- ids2indices(gene_sets, rownames(expr_matrix), remove.empty = TRUE)
idx <- idx[lengths(idx) >= 10 & lengths(idx) <= 500]                  # same size filters as GSEA

cam <- camera(expr_matrix, idx, design, contrast = 2, inter.gene.cor = NA)   # NA = estimate the correlation per set
head(cam[cam$FDR < 0.05, ])          # Direction Up/Down, PValue, FDR
set.seed(123)                        # fry is deterministic; set.seed matters only for roast/mroast
fr <- fry(expr_matrix, idx, design, contrast = 2)
```

Bare `camera(...)` keeps the preset `inter.gene.cor = 0.01` and called 7/20 of the same co-regulated nulls; never omit the argument. CAMERA and fry take the matrix, not a ranked vector, and do not report NES or a leading edge. To also get leading-edge genes, run preranked GSEA on limma `t` as a second view and treat its FDR as screening.

## Per-Sample Scores: ssGSEA and GSVA

**Goal:** Convert an expression matrix into a gene-set-by-sample activity matrix to feed clustering, survival, or a classifier - GSVA itself reports no per-pathway p-value.

**Approach:** GSVA >= 1.50 uses a PARAMETER-OBJECT API: build `gsvaParam(...)` or `ssgseaParam(...)` and pass it to `gsva()`. The old `gsva(expr, gset.idx.list, method=)` signature is defunct.

```r
library(GSVA)
# gene_sets: named list of symbol vectors matching rownames(expr_matrix); minSize/maxSize default to 1/Inf, so set them
gsva_scores  <- gsva(gsvaParam(expr_matrix, gene_sets, kcdf = 'Gaussian', minSize = 10, maxSize = 500))  # unsupervised per-sample set scores
ssgsea_scores <- gsva(ssgseaParam(expr_matrix, gene_sets, minSize = 10, maxSize = 500))                  # ssGSEA via the same dispatch

# per-set contrast on the scores: limma on the score matrix (design as in "Matrix + Design")
tt <- topTable(eBayes(lmFit(gsva_scores, design)), coef = 2, number = Inf)
```

`kcdf` exists only on `gsvaParam` (`ssgseaParam` is rank-based and has none). Its default `'auto'` picks the kernel silently, and it changes GSVA scores between input types: `'Gaussian'` for continuous log-scale data (logCPM, microarray, VST), `'Poisson'` for integer counts. Set it explicitly and log it with the other reporting items above. Checked on GSVA 2.0.7 (signature read from `args(gsvaParam)`; the old `method=` form is defunct).

Use CAMERA/fry or preranked GSEA for a set-level CONTRAST; use ssGSEA/GSVA for a per-sample activity matrix for downstream modeling. limma on GSVA scores also works as a contrast (checked: planted sets FDR 3e-21 and 4e-21, 1 of 20 co-regulated nulls called), but the p-values then describe the scores, not a gene-set null.

## Per-Method Failure Modes

### Ranking by raw p-value
**Trigger:** `gene_list <- -log10(de$pvalue)` with no `sign()`. **Mechanism:** the magnitude is symmetric, so up- and down-regulated genes both land at the top. **Symptom:** NES signs are meaningless; "enriched" sets mix directions. **Fix:** rank by `sign(log2FC) * -log10(pmax(p, 1e-30))`, or use DESeq2 `stat` / limma `t`.

### Preranked p-values treated as correlation-honest
**Trigger:** reporting FDR 0.001 from gseGO/fgsea on a co-regulated set. **Mechanism:** gene permutation assumes gene independence; correlated sets inflate the set-statistic variance, so p is too small. **Symptom:** "significant" pathways that are co-expression and do not replicate. **Fix:** state the permutation type; when the matrix and design exist, make CAMERA with the correlation estimated the default set-level test, `limma::camera(..., inter.gene.cor = NA)` (Wu & Smyth 2012; code in "Matrix + Design"; on 20 co-regulated null sets preranked GSEA called 9, CAMERA `NA` 0) - the default preset `inter.gene.cor = 0.01` assumes low correlation and gives no protection against exactly this failure mode (audit-verified: FDR 7.2e-08 on a set with true inter-gene correlation ~0.31 under the default preset, vs 0.80 once estimated).

### Unsorted or duplicated geneList
**Trigger:** an un-sorted vector, or duplicate gene names after ID conversion. **Mechanism:** clusterProfiler validates both before running, not after. **Symptom:** a hard error, not a silent miscalculation - unsorted: `geneList should be a decreasing sorted vector...`; duplicated names: `Duplicate values in names(stats) not allowed`. **Fix:** `sort(gl[!duplicated(names(gl))], decreasing=TRUE)`; prefer a continuous metric (Wald stat / moderated t rarely tie).

### Bare log2FC ranking
**Trigger:** ranking by `log2FoldChange` from raw counts. **Mechanism:** a gene with 2 vs 8 counts shows a huge unstable LFC. **Symptom:** the leading edge is one or two low-count outliers, not a coordinated shift. **Fix:** rank by `stat`/`t`; if LFC is unavoidable use apeglm/ashr-shrunk LFC (never `type='normal'`).

### No set.seed
**Trigger:** running gseGO without fixing the seed. **Mechanism:** the multilevel Monte Carlo is stochastic. **Symptom:** p-values and the significant-set list drift across identical reruns. **Fix:** `set.seed(123)` AND `seed=TRUE` in the call.

### Tiny leading edge believed
**Trigger:** trusting a high |NES| without inspecting `core_enrichment`. **Mechanism:** a 1-2 gene leading edge is outlier-driven, not a pathway shift; large sets reach high |NES| by chance. **Symptom:** an unreplicated headline pathway. **Fix:** FDR first, then leading-edge size/concentration, then NES for prioritization.

### Stale nPerm / FDR<0.25 lore
**Trigger:** copying `nPerm=10000` and "FDR < 0.25" from a Broad-desktop tutorial. **Mechanism:** `nPerm` is silently ACCEPTED, not rejected - clusterProfiler emits two warnings ("We do not recommend using nPerm...", "You are trying to run fgseaSimple...") and falls back from `fgseaMultilevel` to the coarser `fgseaSimple` engine; clusterProfiler `p.adjust` is BH, not the Broad empirical-null FDR that 0.25 was calibrated for. **Symptom:** no error - a complete-looking `gseaResult` computed on the coarser engine, plus a mis-transplanted threshold (audit-verified on clusterProfiler 4.14.6: `nPerm=10000` runs to completion with only warnings). **Fix:** drop `nPerm` entirely and govern tiny-p with `eps`; if the call isn't fully under your control, guard it: `if ('nPerm' %in% names(gse@params)) stop("nPerm forced a fgseaSimple fallback - remove it")`.

### Stale GSVA call
**Trigger:** `gsva(expr, gene_sets, method='ssgsea')`. **Mechanism:** GSVA >= 1.50 dispatches on a parameter object's class. **Symptom:** the old signature errors. **Fix:** `gsva(gsvaParam(expr, gene_sets))` / `ssgseaParam(...)` (a `method=` error means this).

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| `exponent = 1` | Subramanian 2005 *PNAS* 102:15545 | weights each hit by |statistic|; p=0 is the unweighted KS that flags middle-clustered sets with no strong genes |
| `minGSSize = 10` | clusterProfiler default | drops tiny sets that overfit on one outlier |
| `maxGSSize = 500` | clusterProfiler default | drops overly broad sets that always 'enrich' |
| `eps = 0` (default 1e-10) | clusterProfiler / fgsea | replaces nPerm: `eps=0` resolves exact tiny p-values; 1e-10 is the default floor |
| `pAdjustMethod = 'BH'` | clusterProfiler default | Benjamini-Hochberg FDR; NOT the Broad empirical-null FDR, so do not reflex to 0.25 (reserve it for exploratory work and say which cutoff was used) |
| `pvalueCutoff = 0.05` | clusterProfiler default | filters on p.adjust by default; defensible BH cutoff |
| ~>=7 samples/group | Broad GSEA docs | minimum for a non-degenerate phenotype-permutation null |
| `set.seed(123)` | reproducibility | any fixed seed; the point is to fix the multilevel Monte Carlo |
| `pmax(p, 1e-30)` clamp, not `1e-300` | audit-verified | at `exponent=1` the clamp is a hit weight: `1e-300` gives an exact-zero-p gene ~75x the weight of the strongest measured gene (>16% of total \|stat\| from 4/14,000 genes) and can erase a real signal; `1e-30` keeps the distortion small (~2%) |
| genome-wide `gseGO` runtime | audit-measured | ~155s for 14,000 genes against full GO BP on a single core in this environment; budget minutes, not seconds |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Error about names / wrong ES | geneList not named or not sorted decreasing | `sort(setNames(v, ids), decreasing=TRUE)` |
| `--> No gene can be mapped` (often followed by `'organism' is not a slot in class "NULL"`) | wrong keyType/OrgDb, or non-ENTREZ IDs | `bitr` to the expected ID type first |
| gseKEGG returns 0 terms | ENSEMBL/SYMBOL passed, wrong organism code, or KEGG API down | convert to kegg-id/ENTREZ; check the `organism` code; retry (live API) |
| Different results each run | no `set.seed`, or live KEGG DB changed | fix the seed; pin the KEGG run date |
| fgsea warning `There are ties in the preranked stats (X% of the list)` | duplicated or low-precision ranking values (rounded stats, integer ranks, many exact zeros) | a high X (checked: 93% on a stat rounded to 1 decimal) means a degenerate metric - re-derive a continuous `stat`/`t`; a few percent on a continuous Wald/t ranking is benign |
| No enriched terms (no error) | unsigned metric (directions cancel), or the `minGSSize`/`pvalueCutoff` filter leaves nothing | confirm the ranking is signed; rerun with `pvalueCutoff = 1` to see all sets and their `p.adjust` |
| `core_enrichment` is NA / all-ID | `setReadable` not applied | `setReadable(gse, OrgDb, keyType='ENTREZID')` |

## References

- Mootha VK, Lindgren CM, Eriksson KF, et al. 2003. PGC-1alpha-responsive genes involved in oxidative phosphorylation are coordinately downregulated in human diabetes. *Nat Genet* 34:267-273.
- Subramanian A, Tamayo P, Mootha VK, et al. 2005. Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles. *PNAS* 102:15545-15550.
- Korotkevich G, Sukhov V, Budin N, et al. Fast gene set enrichment analysis. *bioRxiv* 060012 (preprint). DOI 10.1101/060012.
- Wu D, Smyth GK. 2012. Camera: a competitive gene set test accounting for inter-gene correlation. *Nucleic Acids Res* 40:e133.
- Wu D, Lim E, Vaillant F, et al. 2010. ROAST: rotation gene set tests for complex microarray experiments. *Bioinformatics* 26:2176-2182.
- Barbie DA, Tamayo P, Boehm JS, et al. 2009. Systematic RNA interference reveals that oncogenic KRAS-driven cancers require TBK1. *Nature* 462:108-112.
- Hanzelmann S, Castelo R, Guinney J. 2013. GSVA: gene set variation analysis for microarray and RNA-seq data. *BMC Bioinformatics* 14:7.
- Reimand J, Isserlin R, Voisin V, et al. 2019. Pathway enrichment analysis and visualization of omics data using g:Profiler, GSEA, Cytoscape and EnrichmentMap. *Nat Protoc* 14:482-517.
- Wijesooriya K, Jadaan SA, Perera KL, et al. 2022. Urgent need for consistent standards in functional enrichment analysis. *PLoS Comput Biol* 18:e1009935.

## Related Skills

- go-enrichment - Gene-list ORA alternative when no ranking exists
- kegg-pathways - KEGG pathway and module enrichment and ID conventions
- reactome-pathways - Reactome curated-pathway enrichment (local DB)
- wikipathways - WikiPathways community-pathway enrichment
- enrichment-visualization - gseaplot2, ridgeplot, and dotplot of GSEA results
- differential-expression/de-results - Source of the ranking statistic and the padj column conventions
- workflows/expression-to-pathways - End-to-end DE-to-enrichment pipeline
