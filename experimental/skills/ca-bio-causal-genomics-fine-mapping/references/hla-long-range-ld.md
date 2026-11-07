# HLA and Long-Range LD

## HLA and Long-Range LD: When to Stop

The HLA region (chr6:28-34 Mb), chromosome 8 inversion (chr8:8-12 Mb), and a handful of other extended LD blocks violate the assumptions of every fine-mapping method.

**Symptoms of irrecoverable LD structure:** Credible sets contain 30+ SNPs at low purity even with L=30; SuSiE-inf credible sets remain wide; `kriging_rss` flags hundreds of SNPs.

**Options:**
- Stratify by classical HLA allele (HIBAG, SNP2HLA imputation) and test allelic series
- Conditional analysis on the lead variant before fine-mapping the residual
- Exclude the region from genome-wide fine-mapping summaries and report separately
- For chr8 inversion: stratify by inversion genotype if known

Document the caveat in any methods section; standard PIPs at HLA are not interpretable as causality estimates.
