# Anticipated Reviewer Pushback

## Anticipated Reviewer Pushback

| Pushback | Standard response |
|----------|-------------------|
| "In-sample vs reference LD?" | In-sample preferred when cohort genotypes available; if reference, report `estimate_s_rss` lambda < 0.05 plus `kriging_rss` outlier count |
| "Credible-set purity?" | `min_abs_corr >= 0.5` (r2 >= 0.25) default; reported per set; relaxed only with explicit rationale for rare-variant fine-mapping |
| "Is L set high enough?" | If returned CS count < L cap: OK (susieR auto-prunes); otherwise raise L. HLA needs L=20-30 |
| "Why not SuSiE-inf?" | Polygenic-shoulder test: count SNPs with marginal -log10(p) > 4 outside the lead credible set; > 50 indicates a polygenic shoulder and SuSiE-inf (Cui 2024) should be used |
| "Why no functional priors?" | PolyFun applied (or manual coding-variant prior used) and reported; if uniform, justify (low-N, mismatched-ancestry baseline-LF) |
| "Credible set has 50 SNPs -- is that fine-mapping?" | Acknowledged as imprecise; reported alongside diagnostics; cross-trait colocalization or functional fine-mapping (PolyFun, MPRA, allelic series) recommended for resolution |
| "Was Neff used for case-control?" | Yes: `Neff = 4/(1/Ncase + 1/Ncontrol)`; report the value used |
| "Allele harmonization?" | Yes: flipped z when GWAS effect allele differs from reference A1; palindromic SNPs at MAF > 0.42 dropped |
