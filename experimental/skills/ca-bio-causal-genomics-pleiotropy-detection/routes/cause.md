# CAUSE
~~~sh
Rscript scripts/cause_fit.R merged.tsv nuisance_snps.txt signal_snps.txt cause.rds
~~~
Supply a genome-wide merged table with snp, beta_hat_1, seb1, beta_hat_2, seb2, A1 and A2, and one SNP ID per line in each list. The nuisance list is sampled from genome-wide variants; the signal list is LD-pruned and selected for exposure. Prepare and validate these lists for the study. The source recommends at least 100 signal variants. Check installed CAUSE and loo interfaces and comparison diagnostics in [CAUSE notes](../references/cause.md).
**Done when:** the saved models can be inspected for fit, shared-factor estimate and causal versus sharing model comparison.
