# Functional priors

```bash
polyfun.py --compute-h2-L2 --no-partitions --output-prefix polyfun_h2 --sumstats gwas_munged.sumstats --ref-ld-chr UKB_baseline_LF/baselineLF2.2.UKB. --w-ld-chr UKB_baseline_LF/weights.UKB.
```

Use ancestry-matched annotation and align SNPVAR to the locus, then pass normalized weights as prior_weights, not prior_variance, as in `references/polyfun-functional-priors.md`. Compare uniform prior results and state prior source.

Done when weighted and uniform results and their assumptions are reported.
