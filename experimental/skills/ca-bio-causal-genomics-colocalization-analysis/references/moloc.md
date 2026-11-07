# moloc

## moloc Multi-Omic Framework (3-5 Traits)

For k traits, moloc tests `2^k - 1` hypotheses. 3 traits -> 15 hypotheses (H_a, H_b, H_c, H_ab, H_ac, H_bc, H_abc, plus "none of the above"); 4 traits -> 31; 5 traits -> 63. The hypothesis H_{all-share} (all k share a single causal variant) is the multi-omic analog of PP.H4.

```r
# moloc 3-trait example; install via remotes::install_github('clagiamba/moloc')
library(moloc)
# Input: list of k dataframes with BETA, SE, N, MAF per SNP and shared SNP IDs
result_moloc <- moloc_test(listData=list(gwas=gwas_df, eqtl=eqtl_df, sqtl=sqtl_df),
                            prior_var=c(0.01, 0.1, 0.5), priors=c(1e-4, 1e-6, 1e-7))
# PPA: posterior over all 15 hypotheses (3-trait case)
# Key column: PPA.abc (all-three-share)
```

moloc is computationally tractable up to k = 5 but explodes beyond; use HyPrColoc for k >= 6.
