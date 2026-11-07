# Quantitative Thresholds
| Quantity | Threshold | Source / Rationale |
|----------|-----------|-------------------|
| Credible set coverage (well-powered GWAS) | 0.95 (default) | Wang 2020 JRSSB; standard convention |
| Credible set coverage (low-N eQTL, GTEx tissue) | 0.9 or 0.8 | At N < 1000, default 0.95 returns very wide CS; document choice in methods |
| Credible set purity (rare-variant fine-mapping) | min_abs_corr >= 0.1 | LD genuinely sparse; relax to retain signal |
| Credible set purity (default common-variant) | min_abs_corr >= 0.5 (r2 >= 0.25) | susieR default; below this the set is LD-confounded |
| Credible set purity (publication-strict) | min_abs_corr >= 0.7 | Stringent claim; rare in practice |
| PIP suggestive | > 0.5 | Convention; "more likely than not causal among set" |
| PIP strong | > 0.9 | Convention; high-confidence single candidate |
| PIP very strong | > 0.95 | Convention; near-certain candidate within credible set |
| L (default cap) | 10 | susieR default; sufficient for most non-HLA loci |
| L (HLA / complex loci) | 20-30 | Empirical; HLA hosts > 10 independent signals for many traits |
| `n` for case-control susie_rss | Neff = 4/(1/Ncase + 1/Ncontrol), NOT Ntotal | Privé F et al 2022 HGG Adv 3:100136; matches the SE scale of logistic-regression sumstats |
| `estimate_s_rss` lambda acceptable | < 0.05 | susieR vignette; > 0.10 indicates serious LD mismatch |
| `kriging_rss` per-SNP flag | abs(`z_std_diff`) > 3 (`|z_obs - z_exp|` elsewhere in this file) | susieR vignette; flag for manual review |
| Locus window (default) | +/- 500 kb from sentinel | Conventional; covers most LD blocks |
| Locus window (conditional-p floor) | Extend until conditional -log10(p) < 4 | Avoids truncating a secondary signal whose conditional evidence leaks into the window edge |
| Locus window (long-range LD) | 5+ Mb or stratify | HLA chr6:25-35Mb, chr8 inversion chr8:8.1-11.9Mb hg38, chr17 H1/H2 inversion |
| FINEMAP `--n-causal-snps` | 5 | Default; raise for HLA |
| FINEMAP `--prob-conv-sss-tol` | 0.001 | Convergence tolerance; rarely needs change (checked on FINEMAP 1.4.2) |
