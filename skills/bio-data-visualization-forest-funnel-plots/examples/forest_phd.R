# Runnable reference: metafor 4.4+, survival, survminer, MendelianRandomization.
# The meta-analysis and Cox data below are labelled synthetic demonstrations, not clinical findings.

suppressPackageStartupMessages({
  library(metafor)
  library(survival)
  library(survminer)
  library(MendelianRandomization)
})
dir.create("plots", showWarnings = FALSE)

# 1. Synthetic study-level log-ORs with visible heterogeneity (k = 12).
studies <- data.frame(
  author = paste("Synthetic study", seq_len(12)), year = 2010:2021,
  log_or = log(c(0.48, 0.62, 0.71, 0.83, 0.95, 1.08, 1.20, 1.42, 1.61, 1.82, 2.05, 2.28)),
  log_or_se = c(0.17, 0.16, 0.20, 0.15, 0.18, 0.16, 0.19, 0.17, 0.21, 0.18, 0.20, 0.22)
)
stopifnot(nrow(studies) >= 10)
res <- rma(yi = log_or, vi = log_or_se^2, data = studies,
           slab = paste(author, year), method = "REML")

meta_label <- function(x) {
  if (x$k < 5) sprintf("RE model (k = %d; HKSJ CI; heterogeneity withheld)", x$k) else
    sprintf("RE model (tau^2 = %.3f; I^2 = %.1f%%; Q p = %s)",
            x$tau2, x$I2, format.pval(x$QEp, digits = 2))
}
pred <- predict(res)
ticks <- log(c(0.25, 0.5, 1, 2, 4))
study_lb <- studies$log_or - qnorm(0.975) * studies$log_or_se
study_ub <- studies$log_or + qnorm(0.975) * studies$log_or_se
limits <- range(c(study_lb, study_ub, res$ci.lb, res$ci.ub,
                  pred$pi.lb, pred$pi.ub, ticks), finite = TRUE)
pdf("plots/forest.pdf", width = 8, height = 6)
forest(res, atransf = exp, at = ticks, alim = limits, refline = 0,
       xlab = "Odds ratio (95% CI)", header = c("Study", "OR [95% CI]"),
       mlab = meta_label(res), addpred = TRUE)
dev.off()
cat(sprintf("Synthetic REML OR %.2f [%.2f, %.2f]; tau^2 %.3f; I^2 %.1f%%; Q p %s\n",
            exp(res$b[1]), exp(res$ci.lb), exp(res$ci.ub), res$tau2, res$I2,
            format.pval(res$QEp, digits = 3)))

# 2. Funnel diagnostics: Egger is guarded at k >= 10; trim-and-fill is sensitivity-only.
pdf("plots/funnel.pdf", width = 7, height = 5)
funnel(res, level = c(90, 95, 99), shade = c("white", "gray55", "gray75"),
       refline = as.numeric(coef(res)[1]), legend = TRUE, xlab = "log(OR)")
dev.off()
if (res$k >= 10) {
  egger <- regtest(res, model = "lm", predictor = "sei")
  cat("Synthetic Egger p =", format.pval(egger$pval, digits = 3), "(asymmetry, not proof of publication bias)\n")
} else message("Egger withheld: k < 10")
res_tf <- trimfill(res)
cat(sprintf("Synthetic trim-and-fill: %d imputed; original OR %.2f; sensitivity OR %.2f\n",
            res_tf$k - res$k, exp(res$b[1]), exp(res_tf$b[1])))

# 3. Synthetic Cox subgroup display. The interaction p tests effect modification.
set.seed(20260924)
n <- 180
clinical_df <- data.frame(
  time = rexp(n, rate = 0.004),
  event = rbinom(n, 1, 0.68),
  treatment = factor(rbinom(n, 1, 0.5), labels = c("control", "treated")),
  subgroup = factor(rep(c("A", "B", "C"), each = n / 3)),
  age = round(rnorm(n, 62, 8))
)
if (any(table(clinical_df$subgroup) < 20) || any(tapply(clinical_df$event, clinical_df$subgroup, sum) < 5)) {
  stop("Synthetic subgroup data unexpectedly sparse")
}
main_fit <- coxph(Surv(time, event) ~ treatment + subgroup + age, data = clinical_df)
int_fit <- coxph(Surv(time, event) ~ treatment * subgroup + age, data = clinical_df)
interaction_p <- anova(main_fit, int_fit, test = "LRT")[2, "Pr(>|Chi|)"]

sub_hr <- do.call(rbind, lapply(levels(clinical_df$subgroup), function(g) {
  d <- droplevels(clinical_df[clinical_df$subgroup == g, ])
  f <- coxph(Surv(time, event) ~ treatment + age, data = d)
  ci <- confint(f)[1, ]
  c(log_hr = unname(coef(f)[1]), lo = unname(ci[1]),
    hi = unname(ci[2]), n = nrow(d))
}))
pdf("plots/cox_subgroup_forest.pdf", width = 7, height = 5)
forest(x = sub_hr[, "log_hr"], ci.lb = sub_hr[, "lo"], ci.ub = sub_hr[, "hi"],
       slab = paste0(levels(clinical_df$subgroup), " (n=", sub_hr[, "n"], ")"),
       atransf = exp, refline = 0, xlab = "Treatment hazard ratio (95% CI)")
title(sprintf("Synthetic treatment-by-subgroup interaction p = %s", format.pval(interaction_p, digits = 3)))
dev.off()
pdf("plots/cox_adjusted_covariates.pdf", width = 9, height = 5)
print(ggforest(main_fit, data = clinical_df, main = "Synthetic adjusted covariate HRs",
               cpositions = c(0.02, 0.22, 0.4), fontsize = 0.7, noDigits = 2))
dev.off()
cat("Synthetic interaction p =", format.pval(interaction_p, digits = 3), "\n")

# 4. MR method comparison on bundled lipid/CHD data. Suppress SNP rows so an outlier cannot hide method diamonds.
mr_dat <- mr_input(bx = ldlc, bxse = ldlcse, by = chdlodds, byse = chdloddsse)
mr_results <- mr_allmethods(mr_dat)
print(mr_results)
mr_egger_result <- mr_egger(mr_dat)
cat(sprintf("MR-Egger intercept %.4f; p = %s\n", mr_egger_result@Intercept,
            format.pval(mr_egger_result@Pleio.pval, digits = 3)))
pdf("plots/mr_method_forest.pdf", width = 8, height = 4)
print(mr_forest(mr_dat, snp_estimates = FALSE,
                methods = c("ivw", "wmedian", "mbe", "egger")))
dev.off()
