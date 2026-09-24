---
name: bio-data-visualization-forest-funnel-plots
description: Build and audit forest plots for study-level effects and funnel plots for meta-analysis asymmetry, including random-effects pooling, small-study safeguards, Cox covariate versus subgroup displays, and Mendelian-randomization method comparisons. Use for meta-analysis, subgroup interaction, funnel, Egger, trim-and-fill, or MR forest-plot requests.
tool_type: r
primary_tool: metafor
license: MIT
author: GPTomics
---

# Forest and Funnel Plots

Use a forest plot to show estimates and confidence intervals; use a funnel plot to inspect small-study asymmetry. These plots describe evidence. They do not establish a treatment recommendation, publication bias, effect modification, or causality on their own.

## Choose the right display

| Question | Use | Do not call it |
|---|---|---|
| Pool study-level OR, RR, HR, or beta estimates | `metafor::rma()` plus `forest()` | a single common effect when heterogeneity remains unexplained |
| Show adjusted covariate HRs from one Cox model | `survminer::ggforest()` | a subgroup-treatment forest |
| Compare treatment effects across subgroups | a `treatment * subgroup` Cox model, stratum-specific estimates, and an interaction test | significant effect modification from visual separation alone |
| Diagnose funnel asymmetry | funnel, contour funnel, and (only at adequate k) Egger | proof of publication bias |
| Compare MR methods | `MendelianRandomization::mr_forest()` | a causal conclusion from one method |

Ratios must use a log x-axis. For a random-effects meta-analysis, report the pooled estimate and CI, tau-squared, Q-test, I-squared only when `k >= 5`, and a prediction interval when it is informative. I-squared is not stable enough to report below five studies. For `k < 3`, do not pool.

## Reference REML forest

`at` controls tick positions and can also constrain the visual range. Set `alim` from the observed intervals so a tick preference does not clip confidence or prediction intervals.

```r
library(metafor)

# studies: author, year, log_or, log_or_se (or analogous log-HR / log-RR columns)
if (nrow(studies) < 3) stop("Do not pool fewer than 3 studies; plot individual estimates only.")

res <- rma(yi = log_or, vi = log_or_se^2, data = studies,
           slab = paste(author, year), method = "REML",
           test = if (nrow(studies) < 5) "knha" else "z")

heterogeneity_label <- function(x) {
  if (x$k < 5) {
    sprintf("RE model (k = %d; HKSJ CI; heterogeneity metrics withheld)", x$k)
  } else {
    sprintf("RE model (tau^2 = %.3f; I^2 = %.1f%%; Q p = %s)",
            x$tau2, x$I2, format.pval(x$QEp, digits = 2))
  }
}

ticks <- log(c(0.25, 0.5, 1, 2, 4))
pred <- predict(res)
study_lb <- studies$log_or - qnorm(0.975) * studies$log_or_se
study_ub <- studies$log_or + qnorm(0.975) * studies$log_or_se
limits <- range(c(study_lb, study_ub, res$ci.lb, res$ci.ub,
                  pred$pi.lb, pred$pi.ub, ticks), finite = TRUE)

forest(res,
       atransf = exp,
       at = ticks,
       alim = limits,
       refline = 0,
       xlab = "Odds Ratio (95% CI)",
       header = c("Study", "OR [95% CI]"),
       mlab = heterogeneity_label(res),
       addpred = TRUE)
```

For `k < 5`, HKSJ is generally preferable to an unadjusted z CI, but a result at `k <= 3` remains highly imprecise; do not describe it as well calibrated or decisive. Investigate substantial heterogeneity (often I-squared above 50%) with pre-specified subgroup analysis or meta-regression before narrating a pooled effect as biologically coherent.

## Funnel, Egger, and trim-and-fill

```r
funnel(res, xlab = "log(OR)", refline = as.numeric(coef(res)[1]))
funnel(res, level = c(90, 95, 99), shade = c("white", "gray55", "gray75"),
       refline = 0, legend = TRUE, xlab = "log(OR)")

if (res$k >= 10) {
  egger <- regtest(res, model = "lm", predictor = "sei")
  message("Egger p = ", format.pval(egger$pval, digits = 3))
} else {
  message("Egger test withheld: k < 10. Inspect the contour funnel descriptively.")
}

res_tf <- trimfill(res)
message(sprintf("Trim-and-fill sensitivity: %d imputed studies; original OR %.2f; adjusted OR %.2f",
                res_tf$k - res$k, exp(res$b[1]), exp(res_tf$b[1])))
```

Treat trim-and-fill as a sensitivity analysis: imputed studies are hypothetical. Funnel asymmetry can also arise from heterogeneity, outcome reporting, chance, or study quality.

## Cox: adjusted covariate forest versus a subgroup forest

`ggforest()` draws adjusted covariate associations from one Cox model. It is not a treatment-by-subgroup display.

```r
library(survival)
library(survminer)

fit_adjusted <- coxph(Surv(time, event) ~ treatment + age + stage, data = clinical_df)
ggforest(fit_adjusted, data = clinical_df, main = "Adjusted covariate hazard ratios",
         cpositions = c(0.02, 0.22, 0.4), fontsize = 0.7,
         refLabel = "Reference", noDigits = 2)
```

For a prespecified subgroup analysis, first reject sparse strata. Then estimate treatment HRs within strata and report the interaction likelihood-ratio p-value; a subgroup-specific CI is not an interaction test.

```r
clinical_df$subgroup <- droplevels(factor(clinical_df$subgroup))
counts <- table(clinical_df$subgroup)
if (any(counts < 20)) stop("Collapse or omit sparse subgroup levels before estimating subgroup HRs.")

main_fit <- coxph(Surv(time, event) ~ treatment + subgroup, data = clinical_df)
int_fit <- coxph(Surv(time, event) ~ treatment * subgroup, data = clinical_df)
interaction_p <- anova(main_fit, int_fit, test = "LRT")[2, "Pr(>|Chi|)"]

by_subgroup <- do.call(rbind, lapply(levels(clinical_df$subgroup), function(g) {
  d <- droplevels(clinical_df[clinical_df$subgroup == g, ])
  if (sum(d$event) < 5) stop("Subgroup ", g, " has fewer than 5 events.")
  f <- coxph(Surv(time, event) ~ treatment, data = d)
  ci <- confint(f)[1, ]
  c(group = g, log_hr = unname(coef(f)[1]), lo = unname(ci[1]),
    hi = unname(ci[2]), n = nrow(d))
}))
by_subgroup <- as.data.frame(by_subgroup, stringsAsFactors = FALSE)
by_subgroup[, c("log_hr", "lo", "hi", "n")] <- lapply(by_subgroup[, c("log_hr", "lo", "hi", "n")], as.numeric)

forest(x = by_subgroup$log_hr, ci.lb = by_subgroup$lo, ci.ub = by_subgroup$hi,
       slab = paste0(by_subgroup$group, " (n=", by_subgroup$n, ")"),
       atransf = exp, refline = 0, xlab = "Treatment hazard ratio (95% CI)")
title(sprintf("Treatment-by-subgroup interaction p = %s", format.pval(interaction_p, digits = 3)))
```

## MR method-comparison forest

Outlying SNP estimates can make method diamonds unreadable. Use `snp_estimates = FALSE` for the method comparison, report SNP-level diagnostics separately, and print method results rather than promising p-values on a plot that does not support them.

```r
library(MendelianRandomization)

mr_dat <- mr_input(bx = bx, bxse = bxse, by = by, byse = byse)
mr_results <- mr_allmethods(mr_dat)
print(mr_results)
egger <- mr_egger(mr_dat)
message(sprintf("MR-Egger intercept = %.4f (p = %s)",
                egger@Intercept, format.pval(egger@Pleio.pval, digits = 3)))
mr_forest(mr_dat, snp_estimates = FALSE,
          methods = c("ivw", "wmedian", "mbe", "egger"))
```

## Optional layouts

`forestplot` is useful when table columns matter; supply `boxsize` so weight is visible. `netmeta` is a frequentist network-meta-analysis package, not a Bayesian one. `cumul()` makes a cumulative meta-analysis in temporal order.

```r
# forestplot::forestplot(labeltext = table_text, mean = estimates,
#   lower = lower_ci, upper = upper_ci, zero = 1, xlog = TRUE,
#   boxsize = sqrt(weights(res)) / max(sqrt(weights(res))))

# net <- netmeta::netmeta(TE, seTE, treat1, treat2, studlab, sm = "OR", random = TRUE)
# netmeta::forest(net, random = TRUE)

# cumulative <- cumul(res, order = studies$year)
# forest(cumulative, atransf = exp, refline = 0)
```

## Failure checks before reporting

- Do not pool `k < 3`; use HKSJ and withhold I-squared below five studies.
- Do not run or interpret Egger below ten studies.
- Keep all CIs and the prediction interval inside `alim`; arrows indicate clipping, not a complete interval.
- Use a log axis for OR, RR, and HR.
- Present original and trim-and-fill estimates together, labeling the latter sensitivity-only.
- Require a formal interaction test for subgroup claims and collapse sparse strata.
- Triangulate MR methods and report MR-Egger's intercept; do not infer causality from plot agreement alone.

## References

- Egger M, et al. 1997. *BMJ* 315:629-634.
- Sterne JAC, et al. 2011. *BMJ* 343:d4002.
- Duval S, Tweedie R. 2000. *Biometrics* 56:455-463.
- Higgins JPT, Thompson SG. 2002. *Stat Med* 21:1539-1558.
- Higgins JPT, Thompson SG, Spiegelhalter DJ. 2009. *JRSS-A* 172:137-159.
- IntHout J, Ioannidis JPA, Borm GF. 2014. *BMC Med Res Methodol* 14:25.
- Borenstein M, et al. 2017. *Res Synth Methods* 8:5-18.
- Bowden J, Davey Smith G, Burgess S. 2015. *Int J Epidemiol* 44:512-525.
