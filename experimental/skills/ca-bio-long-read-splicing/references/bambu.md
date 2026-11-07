## Bambu for Annotation-Aware Discovery + Quantification

**Goal:** Joint discovery and quantification with statistical filtering of novel isoforms.

**Approach:** R Bioconductor package; takes BAM + reference annotation + genome; outputs ranged SE objects of known + novel transcripts.

```r
library(bambu)

bam_files <- c('sample1.bam', 'sample2.bam', 'sample3.bam')
genome <- 'reference.fa'
gtf <- 'gencode.v45.annotation.gtf'

bambuAnnotations <- prepareAnnotations(gtf)

se <- bambu(
    reads = bam_files,
    annotations = bambuAnnotations,
    genome = genome,
    NDR = 0.1,
    ncore = 8   # BiocParallel workers; on Windows R use ncore = 1
)

writeBambuOutput(se, path = 'bambu_output/')

tx_counts <- as.data.frame(assays(se)$counts)
gene_counts <- transcriptToGeneExpression(se)
```

Bambu (Chen 2023 *Nat Methods* 20:1187-1195) uses **NDR** (Novel Discovery Rate) as a single, calibrated parameter replacing per-sample heuristics:

| NDR | Interpretation |
|-----|----------------|
| 0.05 | Stringent; few novel transcripts; highest precision |
| 0.1 | Balanced (default) |
| 0.2-0.3 | Permissive; more novel discoveries; recall over precision |

Excellent for combined discovery + quantification when statistical filtering matters. With several BAMs on Windows R, `ncore = 8` stopped with "BiocParallel errors ... could not find function seqlengths" and `ncore = 1` ran (single-BAM runs worked with either; multicore was not tried on Linux). NDR needs enough read classes: Bambu warns "NDR approximated" below ~50 read classes, and on a 3.4k-read simulated set NDR 0.05 and 0.1 reported no novel transcript (0.3: one; 1.0: two). Use it on genome-scale data, not a few genes.
