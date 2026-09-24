# PolyFun and Manual Functional Priors

## Functional Priors with PolyFun

**Goal:** Use genome-wide stratified LDSC heritability to weight per-SNP causal priors, sharpening PIPs at coding, conserved, and regulatory variants.

**Approach:** Run PolyFun once genome-wide to estimate per-SNP h2 from the baseline-LF annotation set; extract per-SNP causal prior; pass to susie_rss as `prior_weights`.

```bash
# Parametric route: L2-regularized S-LDSC writes per-SNP priors directly (--no-partitions)
polyfun.py --compute-h2-L2 --no-partitions \
    --output-prefix polyfun_h2 \
    --sumstats gwas_munged.sumstats \
    --ref-ld-chr UKB_baseline_LF/baselineLF2.2.UKB. \
    --w-ld-chr UKB_baseline_LF/weights.UKB.
# Per-SNP priors written to polyfun_h2.<CHR>.snpvar_ridge_constrained.gz

# Non-parametric route (finer, optional): drop --no-partitions above, then add an
# intermediate LD-score step before re-estimating binned per-SNP h2:
#   polyfun.py --compute-ldscores --output-prefix polyfun_h2 ...
#   polyfun.py --compute-h2-bins --output-prefix polyfun_h2 --sumstats gwas_munged.sumstats --w-ld-chr UKB_baseline_LF/weights.UKB.
```

```r
library(susieR)
priors <- read.table('polyfun_h2.6.snpvar_ridge_constrained.gz', header = TRUE)
priors <- priors[match(gwas_df$SNP, priors$SNP), ]
prior_w <- priors$SNPVAR / sum(priors$SNPVAR, na.rm = TRUE)

fit <- susie_rss(z = z_scores, R = ld_matrix, n = N, L = 10,
                 prior_weights = prior_w)
```

UKB baseline-LF priors are pre-computed EUR-only at `data.broadinstitute.org/alkesgroup/UKBB_LD/` for hg19 and hg38. For EAS, AFR, or SAS GWAS, the EUR weights are NOT valid: functional-prior fine-mapping in a non-EUR ancestry requires baseline-LF annotations matched to that ancestry. For ancestries lacking matched baseline-LF (admixed, under-represented), accept reduced power and run uniform-prior susie_rss; applying EUR weights to non-EUR sumstats produces miscalibrated PIPs that look sharper than reality.

### Manual Coding-Variant Priors Without PolyFun

For postdocs without PolyFun infrastructure or with single-locus inputs, manual annotation-based priors are a reasonable approximation (Hutchinson 2020 Hum Mol Genet 29:R81). As a stated convention, coding variants get ~10x uniform weight; broadly conserved variants ~5x (binned by CADD-PHRED quantile).

```r
build_manual_priors <- function(vep_df, cadd) {
    w <- rep(1, nrow(vep_df))
    w[vep_df$Consequence %in% c('missense_variant', 'stop_gained', 'splice_donor_variant',
                                'splice_acceptor_variant', 'frameshift_variant')] <- 10
    w[cadd >= quantile(cadd, 0.95, na.rm = TRUE)] <- pmax(w[cadd >= quantile(cadd, 0.95, na.rm = TRUE)], 5)
    w / sum(w)
}
fit <- susie_rss(z = z_scores, R = ld_matrix, n = Neff, L = 10, prior_weights = build_manual_priors(vep, cadd))
```

Report the prior construction explicitly; reviewers will ask whether the prior was tuned post hoc.
