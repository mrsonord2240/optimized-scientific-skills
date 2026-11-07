# Variant + Outlier Integration


**Goal:** Connect a candidate splice-altering DNA variant to an RNA-level outlier in the same sample, as a research prioritisation.

**Approach:** Cross-reference SpliceAI hits with FRASER2 outliers in the same sample. `delta_max` is the largest of the four SpliceAI delta scores (DS_AG, DS_AL, DS_DG, DS_DL).

```bash
bcftools query -f '%CHROM\t%POS\t%INFO/SpliceAI\n' spliceai_annotated.vcf > spliceai_raw.tsv
```

```r
library(dplyr)

raw <- read.table('spliceai_raw.tsv', sep = '\t', quote = '', col.names = c('chrom', 'pos', 'spliceai'),
                  stringsAsFactors = FALSE)
raw <- raw[raw$spliceai != '.', ]
# INFO = ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL, one comma-separated entry per ALT allele;
# unscored alleles are '.', which max() ignores
raw$delta_max <- vapply(strsplit(raw$spliceai, ',', fixed = TRUE), function(e)
    suppressWarnings(max(0, as.numeric(unlist(lapply(strsplit(e, '|', fixed = TRUE), `[`, 3:6))), na.rm = TRUE)),
    numeric(1))

fraser_hits <- read.table('fraser_results.tsv', header = TRUE, sep = '\t')
# 'chr21' vs '21' otherwise joins to zero rows without a warning
strip_chr <- function(x) sub('^chr', '', x)
variants <- mutate(raw, chrom = strip_chr(chrom))
fraser_hits <- mutate(fraser_hits, seqnames = strip_chr(seqnames))
stopifnot(any(variants$chrom %in% fraser_hits$seqnames))

confirmed <- variants %>%
    filter(delta_max >= 0.2) %>%
    inner_join(
        fraser_hits %>% filter(sampleID == 'PATIENT_001', padjust < 0.05),
        by = c('chrom' = 'seqnames'),
        relationship = 'many-to-many'
    ) %>%
    filter(abs(pos - start) < 1000 | abs(pos - end) < 1000)
```

A SpliceAI hit with a concordant FRASER2 outlier in the same sample is a candidate for functional follow-up. Whether it counts as functional evidence, and at what strength, is a clinical-laboratory decision (ClinGen SVI splicing recommendations, Walker 2023) outside this Skill. The 1 kb window is a loose heuristic.
