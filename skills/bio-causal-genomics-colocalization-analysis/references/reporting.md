# Reporting and reviewer pushback

## Anticipated Reviewer Pushback

| Pushback | Standard response |
|----------|-------------------|
| "Sensitivity to p12 prior?" | `coloc::sensitivity()` reported; PP.H4 robust across 1e-7 to 1e-5 grid |
| "Why not coloc.susie? Multi-causal possible?" | coloc.abf single-causal assumption stated; if PP.H3 dominant or GCTA-COJO identifies >= 2 independent signals, coloc.susie / SuSiE-based run; reported |
| "LD reference matched?" | In-sample preferred; if reference panel used, `estimate_s_rss(z, R, N)` lambda < 0.05; `kriging_rss` diagnostic clean |
| "PP.H4 = 0.6 is colocalization?" | No -- bands stated: 0.5-0.7 suggestive; >= 0.7 triangulation tier; >= 0.8 standard publication; >= 0.95 industry/clinical |
| "MHC region included?" | chr6:25-35 Mb excluded; HLA-coloc (Butler-Laporte 2024) for classical-allele-level coloc |
| "Ancestry mismatch?" | LD reference ancestry-matched to GWAS; for cross-ancestry use coloc_SuSiEx |
| "Sentinel SNP swap?" | Re-centered window on each trait's lead, joint top, eQTL top; PP.H4 stable within 0.1 |

## Reviewer-Grade Reporting Template

For each colocalization claim, the report should include:

1. **Method and version** (e.g. coloc 5.2.3 coloc.abf, or coloc.susie with SuSiE L=10).
2. **Window definition** (e.g. +/- 500 kb around the GWAS lead rs12345 at chr6:30450000, hg38), and lead-SNP-swap sensitivity (PP.H4 at GWAS lead vs eQTL lead vs joint top).
3. **Priors** p1, p2, p12 used; **sensitivity** plot from `coloc::sensitivity()` and the p12 range over which PP.H4 stays above threshold.
4. **All five posteriors** PP.H0 through PP.H4 (not PP.H4 alone).
5. **Threshold band** the result clears (>= 0.7 triangulation tier, >= 0.75 Open Targets screening, >= 0.80 published, >= 0.90 stringent, >= 0.95 clinical).
6. **LD reference** ancestry, source (1000G phase 3 EUR / in-sample / UKBB), and lambda from `estimate_s_rss` if coloc.susie.
7. **Reference QTL panel** version (e.g. GTEx v8 MASHR-EUR, PredictDB release 2022-01; eQTLGen 2019).
8. **Conditional analysis** GCTA-COJO results if multi-causal; per-credible-set PP if coloc.susie.
9. **Failure-mode caveats** explicitly addressed: MHC excluded, chr 8 inversion excluded, ancestry-matched LD, sdY/s correctly specified, palindromic SNPs handled.
10. **Methods-section H0-H4 prose** describing what each hypothesis means (see usage-guide.md).
