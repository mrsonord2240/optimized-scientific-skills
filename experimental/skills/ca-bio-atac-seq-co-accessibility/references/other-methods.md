# Other-Methods



| Tool | Method | Input | Output | Strength | Fails when |
|------|--------|-------|--------|----------|------------|
| Cicero (Pliner 2018) | Graphical lasso on aggregated cell metacells | scATAC peak-cell matrix + cell trajectory | Peak-pair connection score (-1 to 1) | Original, well-validated; integrates with Monocle3 | Slow on large datasets (">50K cells" is an unverified rule of thumb); sensitive to alpha tuning |
| ArchR getCoAccessibility | Cicero-based; uses ArchR's metacell aggregation | ArchR project | Same as Cicero | Built-in to ArchR pipeline; faster on large datasets | Tied to ArchR; same biology as Cicero |
| SCENIC+ (Bravo 2023) | Multi-step: co-accessibility + motif scoring + RNA correlation | Multiome (ATAC + RNA) or paired | TF-driven enhancer-gene networks | Most comprehensive; multi-modal | Multiome data required; computationally heavy |
| LinkPeaks (Signac) | Pearson correlation of accessibility with paired gene expression | Multiome | Peak-gene linkage score | Direct enhancer-gene from RNA correlation | Multiome-only; not pure ATAC |
| GeneHancer / FANTOM5 / EpiMap | Bulk-derived enhancer-gene reference | None (database lookup) | Pre-computed enhancer-gene pairs | Comprehensive; published references | Cell-type-agnostic; may not match the biology of interest |

Methodology evolves; verify against Pliner 2018 (Cicero), Bravo 2023 (SCENIC+), Nasser 2021 (ABC model alternative for enhancer-gene), and current Hi-C concordance benchmarks.



```r
library(ArchR)
proj <- loadArchRProject('ArchR_out')
proj <- addCoAccessibility(proj, reducedDims='IterativeLSI',
                          k=100, knnIteration=500,
                          maxDist=250000)               # 250 kb cis (wider than the 100 kb default)
co_acc <- getCoAccessibility(proj, corCutOff=0.5,       # Default 0.5 in ArchR; lower for more (calibrate against Hi-C/HiChIP)
                             returnLoops=TRUE)           # TRUE (default) -> SimpleList, loops GRanges in element 1; FALSE -> DataFrame of peak-pair correlations
loops <- co_acc[[1]]                                    # GRanges of loops; ArchR's own help calls the return a GRanges, it is a SimpleList
```

With `returnLoops=TRUE` (the default) ArchR returns a SimpleList whose first element is the loops GRanges (compatible with `GenomicInteractions` for direct overlap with Hi-C loops); `returnLoops=FALSE` instead returns a DataFrame of peak-pair correlations.



SCENIC+ (1.0a2) runs as a Snakemake pipeline (CLI), not a single monolithic Python call; this route was reviewed statically only (not installable in the audit environment). Prepare the inputs first (a pycisTopic cisTopic object, motif-enrichment results, and paired RNA AnnData), then scaffold and run the workflow:

```bash
# Scaffold the pipeline, then edit its config.yaml to point at the cisTopic object,
# motif-enrichment results, and GEX AnnData
scenicplus init_snakemake --out_dir scplus_pipeline/
snakemake --cores 16 --snakefile scplus_pipeline/Snakemake/workflow/Snakefile
# eRegulons (TF + target genes + linked enhancers) are written to the output MuData (scplusmdata.h5mu)
```

SCENIC+ is significantly more complex than Cicero; budget 1-2 days for setup. The benefit is that outputs are TF -> enhancer -> gene triples, not just peak-peak co-accessibility.
