# PWCoCo

## PWCoCo (Conditional Pairwise Coloc)

PWCoCo (Robinson 2022) wraps GCTA-COJO conditional analysis around coloc.abf. For a locus with `k1` independent trait-1 signals and `k2` independent trait-2 signals, PWCoCo runs `k1 * k2` pairwise coloc.abf tests after conditioning each summary statistic on the other independent signals.

**When to use:** When GCTA-COJO has identified >= 2 independent signals in at least one trait and individual-level reference genotypes are available. Particularly suited to bulk eQTL with secondary cis signals.

**Inputs:** Per-trait summary stats (SNP, A1, A2, freq, beta, se, p, N) + plink bfile reference. **Output:** One coloc.abf result per (conditional signal 1, conditional signal 2) pair. Interpret each row as an independent single-signal coloc.

**Caveats:** PWCoCo requires individual-level reference (plink bfile); cannot run on summary stats alone. Collinearity threshold in COJO (default `--cojo-collinear 0.9`) controls how aggressively independent signals are split; lower values fragment, higher values merge.

```bash
# Step 1: GCTA-COJO identifies independent signals at the locus
gcta64 --bfile 1KG_EUR --chr 6 --extract locus.snplist \
       --cojo-file gwas.ma --cojo-slct --out gwas_cojo

# gwas_cojo.jma.cojo lists independent signals (per --cojo-p 5e-8 default)

# Step 2: PWCoCo runs pairwise conditional coloc.abf per signal pair
pwcoco --bfile 1KG_EUR --sum_stats1 gwas.txt --sum_stats2 eqtl.txt \
       --p_cutoff1 5e-8 --p_cutoff2 5e-5 \
       --chr 6 \
       --out pwcoco_result
```

Output: one coloc.abf result per (conditional signal in trait 1, conditional signal in trait 2) pair. Interpret each row as a separate single-signal coloc test. If COJO finds 2 GWAS + 1 eQTL signals, expect 2 result rows.
