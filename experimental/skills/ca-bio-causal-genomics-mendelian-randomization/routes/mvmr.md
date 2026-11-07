# Multivariable MR

```bash
Rscript scripts/mvmr_conditional_f.R --dat mvmr_input.tsv --exposures x1,x2 --gencov 0
```

Input has SNP, beta.<exposure>, se.<exposure>, beta.y and se.y. Zero genetic covariance is a source assumption for non-overlapping exposures; the wrapper accepts only numeric scalar gencov and does not implement an overlap covariance matrix. Stop and adapt the analysis when overlap needs a matrix. Report conditional F for each exposure and Q_A; below 10 is weak and the script stops below 1. See `references/mvmr-conditional-f.md`.

Done when estimates, conditional strengths and overlap assumptions are reported.
