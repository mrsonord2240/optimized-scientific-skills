# LHC-MR and Choosing CAUSE vs LHC-MR

## LHC-MR Workflow

**Goal:** Jointly estimate forward causal effect, reverse causal effect, and the heritable-confounder contribution from genome-wide sumstats (not just significant SNPs).

```r
library(lhcMR)
merged <- merge_sumstats(list(df_x, df_y), c('X', 'Y'), LD.filepath='ldsc/LDscores.txt', rho.filepath='ldsc/LDrho.txt')
sp_list <- calculate_SP(merged, trait.names=c('X', 'Y'), run_ldsc=TRUE, run_MR=TRUE, hm3='ldsc/w_hm3.snplist', ld='ldsc/eur_w_ld_chr/', nStep=2, SP_single=3, SP_pair=50)
res_lhc <- lhc_mr(sp_list, trait.names=c('X', 'Y'), paral_method='lapply', nBlock=200, nCores=4)
```

LHC-MR is computationally heavy (hours on full sumstats) but among the most rigorous CHP-aware estimators when both GWAS are well-powered. Output includes axx, ayy, hxy (confounder effect on each trait), and bidirectional alpha_xy, alpha_yx.

**Choosing CAUSE vs LHC-MR (Darrous 2021):**

| Condition | Preferred method |
|-----------|------------------|
| `>= 100` genome-wide significant SNPs after pruning | CAUSE (Bayesian; CHP-explicit; mature posterior diagnostics) |
| Polygenic exposure with few significant loci | LHC-MR (uses genome-wide signal, not just significant SNPs) |
| Severe sample overlap between exposure and outcome GWAS | LHC-MR (jointly models overlap); CAUSE's rho correction is exposed to misspecification at high overlap |
| Bidirectionality of central interest | LHC-MR (jointly estimates alpha_xy and alpha_yx); CAUSE only models forward |
| Limited compute / quick turnaround | CAUSE (minutes to hours); LHC-MR may be > 24h on full sumstats |

When both apply, report both with the agreement / disagreement explicit in the discussion.
