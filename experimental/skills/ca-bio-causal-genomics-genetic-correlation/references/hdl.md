# HDL: Lower-Variance rg for Independent Samples
**Goal:** Estimate global rg with ~60% lower variance than LDSC when the two GWAS draw from non-overlapping samples.

**Approach:** Format each GWAS as an HDL data frame; supply the UKB-array SVD eigen reference path; pass `N0` (overlapping sample count, 0 for independent).

```r
# remotes::install_github('zhenin/HDL/HDL')
library(HDL)

gwas1 <- data.frame(
    SNP = trait1$rsid,
    A1 = trait1$effect_allele,
    A2 = trait1$other_allele,
    N = trait1$N,
    Z = trait1$beta / trait1$se,
    b = trait1$beta,
    se = trait1$se
)

gwas2 <- data.frame(
    SNP = trait2$rsid,
    A1 = trait2$effect_allele,
    A2 = trait2$other_allele,
    N = trait2$N,
    Z = trait2$beta / trait2$se,
    b = trait2$beta,
    se = trait2$se
)

res <- HDL.rg(
    gwas1.df = gwas1,
    gwas2.df = gwas2,
    LD.path = 'UKB_array_SVD_eigen90_extraction',
    N0 = 0,  # number of overlapping individuals; 0 for independent cohorts. HDL corrects for overlap via N0 and is robust to its misspecification
    output.file = 'hdl_rg.txt'
)

print(res$rg)
print(res$rg.se)
print(res$P)
```

Pre-download the UKB SVD reference (`HDL_documentation.html` -> "How to obtain LD reference panel" link); `eigen90` is the UKB-array (genotyped-SNP) panel used here, while `eigen99` exists only for the imputed-variant panel -- match the panel to the SNP coverage of the input GWAS, not to a speed/precision setting. Do NOT run HDL when sample overlap is unknown or non-trivial; the wrapper does not warn.
