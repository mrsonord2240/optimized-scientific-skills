---
name: bio-causal-genomics-colocalization-analysis
description: Test whether a GWAS signal and a molecular QTL (eQTL/sQTL/pQTL/mQTL) share the same causal variant -- e.g. "does my GWAS lead SNP colocalize with this eQTL, or is it just LD" -- using Bayesian colocalization (coloc.abf, coloc.susie, HyPrColoc, moloc, eCAVIAR, SMR/HEIDI, PWCoCo, SharePro). Use when integrating GWAS with eQTL/sQTL/pQTL/mQTL, distinguishing shared causal variants from LD-driven coincidence, handling allelic heterogeneity, choosing between single-causal vs multi-causal methods, picking PP.H4 thresholds, running sensitivity over p12, or harmonising summary statistics for colocalization.
tool_type: r
primary_tool: coloc
license: MIT
author: GPTomics
---

## Version Compatibility

Reference examples tested with: coloc 5.2.3+, susieR 0.12.35+, hyprcoloc 0.0.2 (GitHub jrs95/hyprcoloc; the GitHub DESCRIPTION version, checked 2026-09-21), SMR 1.3.1+ (CLI, cnsgenomics.com), eCAVIAR 2.2+ (compiled from caviar/eCAVIAR repo), PWCoCo 1.0+ (jwr-git/pwcoco), moloc 0.1+ (clagiamba/moloc), SharePro_coloc 7.0+ (zhwm/SharePro_coloc), R >= 4.1.

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('coloc')`; check `?coloc.abf`, `?coloc.susie`, `?runsusie`
- CLI: `smr --version`, `pwcoco --help`, `sharepro_coloc.py --help`

If code throws AttributeError, NULL list elements, or `Error in coloc.abf: dataset must have...`, introspect the installed package signature and adapt the example rather than retrying.

# Colocalization Analysis

**"Test whether my GWAS signal and an eQTL share the same causal variant"** -> Compute Bayesian posterior probabilities over five hypotheses (H0 neither, H1 trait-1-only, H2 trait-2-only, H3 distinct causal variants, H4 shared causal variant) to discriminate true causal overlap from LD-driven coincidence, then run sensitivity analysis over the p12 prior.

- R (single-causal, fastest): `coloc::coloc.abf(dataset1, dataset2, p12=5e-6)` -> `coloc::sensitivity(res, 'H4 > 0.75')`
- R (multi-causal, needs LD): `runsusie(d1)` -> `runsusie(d2)` -> `coloc.susie(s1, s2)` -> per-credible-set PP
- R (many traits, single-causal cluster): `hyprcoloc::hyprcoloc(effect.est = betas_mat, effect.se = ses_mat, trait.names = ..., snp.id = ...)` -> trait clusters
- CLI (causality vs linkage): `smr --bfile ref --gwas-summary g.ma --beqtl-summary eqtl.besd --out smr` -> SMR p + HEIDI p
- CLI (allelic heterogeneity): eCAVIAR `eCAVIAR -l ld1 -l ld2 -z z1 -z z2 -o out -c 2` -> CLPP per SNP
- CLI (conditional): PWCoCo conditions on each independent signal via GCTA-COJO then runs pairwise coloc.abf

## Algorithmic Taxonomy

| Method | Model | Inputs | Output | Strength | Fails when |
|--------|-------|--------|--------|----------|------------|
| coloc.abf (Giambartolomei 2014) | Single causal variant per locus; Bayesian ABF | beta+varbeta or p+MAF; sample sizes; type/s/sdY | PP.H0-H4 | Fast (~1s/locus), no LD required, mature, widely-cited | 2+ causal variants in moderate LD -> PP.H3 inflates spuriously; assumes a single causal per trait |
| coloc.susie (Wallace 2021) | Multi-causal via SuSiE; per-credible-set pairwise coloc | Summary stats + ancestry-matched LD matrix | PP.H4 per (CS1, CS2) pair | Handles allelic heterogeneity; principled CS framework | Sensitive to LD-mismatch; sample-size-LD mismatch -> spurious credible sets; needs in-sample or matched LD |
| SMR + HEIDI (Zhu 2016) | Tests pleiotropy (one variant -> both traits) vs linkage (two variants in LD) | GWAS .ma; eQTL .besd; LD reference (plink bfile) | SMR p (significance) + HEIDI p (null = shared causal) | Distinguishes shared-causal from linkage at a top SNP; standard for eQTLGen / GTEx integration | Fails to discriminate when LD between causal SNPs > 0.7 (HEIDI loses power); HEIDI requires >= 10 SNPs near top |
| eCAVIAR / CLPP (Hormozdiari 2016) | Fine-mapping-aware; computes Colocalization Posterior Probability per SNP | Z-scores; LD matrices per trait | CLPP per SNP; per-locus sum | Handles allelic heterogeneity natively; per-SNP resolution | Computationally heavy at -c > 3 causal variants; CLPP thresholds debated (0.01 vs 0.1) |
| PWCoCo (Robinson 2022) | Pairwise conditional via GCTA-COJO conditioning | Summary stats + individual-level LD bfile | Per-conditional-signal coloc.abf results | Cleanly handles AH at top GWAS hit + secondary signals | Needs individual-level reference; sensitive to COJO collinearity threshold |
| moloc (Giambartolomei 2018) | Multi-trait extension of coloc.abf (3-5 traits) | Per-trait summary stats | 15 (3-trait) / 31 (4-trait) / 63 (5-trait) hypothesis PPs | First principled multi-omic coloc | Hypothesis count = 2^k - 1 explodes; >= 6 traits computationally infeasible; minimally updated since 2019 |
| HyPrColoc (Foley 2021) | Many-trait cluster-based; iterative branch-and-bound under single-causal | Beta + SE matrices SNPs x traits | Trait clusters sharing a causal variant | Scales to 50+ traits; identifies cluster substructure | Inherits single-causal assumption from coloc.abf; clusters can fragment under AH |
| SharePro_coloc (Zhang 2024) | Variational effect-group joint model | Beta + SE; LD per ancestry | Effect-group level PP | Handles multiple causal signals jointly; faster than coloc.susie at scale | Newer (2024); benchmarks evolving; trickier installation |
| Pullin & Wallace 2025 variant-specific priors | Function-aware p12 (e.g. up-weight coding/promoter SNPs) | Same as coloc.abf + per-SNP prior weights | PP.H0-H4 with non-uniform prior | Improves discovery when functional annotation is informative | Annotation choice is a methodological lever; report sensitivity |

Methodology evolves; verify the current Open Targets Genetics, eQTL Catalogue, and FinnGen colocalization pipelines before locking parameters. Open Targets uses coloc.abf at PP.H4 >= 0.75 with p12 = 1e-5; FinnGen uses coloc.susie at PP.H4 >= 0.8 with in-sample LD.

## Decision Tree by Scenario

| Scenario | Recommended method | Why |
|----------|---------------------|-----|
| GWAS + single-tissue eQTL, top GWAS variant looks single-signal | coloc.abf + sensitivity() | Fast, no LD needed, well-validated; single-causal assumption typically holds at clean loci. Pipeline below; harmonise first: `references/allele-harmonisation.md` |
| GWAS + eQTL, conditional analysis shows 2+ independent signals | coloc.susie OR PWCoCo | Multi-causal handling; coloc.susie if summary-stats LD available, PWCoCo if individual-level reference accessible. See `references/coloc-susie.md`, `references/pwcoco.md` |
| GWAS + multi-tissue eQTL (e.g. all 49 GTEx tissues) | coloc.abf per tissue + HyPrColoc across tissues | Per-tissue PP.H4 gives tissue-specific causality; HyPrColoc identifies tissue clusters sharing the variant. See `references/hyprcoloc.md` |
| GWAS + eQTL + sQTL + mQTL (3-5 omics) | moloc (k <= 5) OR HyPrColoc | moloc gives explicit hypothesis posterior; HyPrColoc scales but loses hypothesis structure. See `references/moloc.md`, `references/hyprcoloc.md` |
| Many GWAS traits at one locus (pleiotropic hub) | HyPrColoc | Designed for many-trait clustering; coloc.abf pairwise scales as k^2. See `references/hyprcoloc.md` |
| Top SNP has only modest GWAS p; is it the same causal as eQTL? | SMR + HEIDI | SMR tests pleiotropy/causality; HEIDI rejects shared-causal -> linkage. See `references/smr-heidi.md` |
| Want per-SNP credibility under allelic heterogeneity | eCAVIAR (CLPP) | Per-SNP CLPP integrates fine-mapping with coloc. See `references/ecaviar-clpp.md` |
| MHC / HLA region (chr6:25-35 Mb) | HLA-coloc (Butler-Laporte 2024) OR exclude MHC | Long-range LD breaks single-causal assumption; standard PP.H4 not interpretable |
| Trans-eQTL / GWAS pair | coloc.abf with p12 lowered to 5e-6 or 1e-6 | Shared causality is biologically rare; default p12=1e-5 over-favours H4 |
| Ancestry-mismatched GWAS vs eQTL | ancestry-matched coloc.susie OR coloc_SuSiEx | LD differs across ancestries; using EUR LD on AFR z-scores produces spurious credible sets. See `references/coloc-susie.md` |
| Very small eQTL (N < 200) | None reliably; flag locus underpowered | All methods report H0/H1/H2 dominance; report PP transparently and gather larger reference (eQTLGen N~31k, GTEx v8) |

## Per-Method Failure Modes

### coloc.abf -- PP.H3 inflation under multiple causal variants

**Trigger:** Locus has 2+ independent causal signals in moderate LD (r2 ~ 0.3-0.6), one shared with the eQTL and one GWAS-only, **and both datasets are powered to see the signals** (eQTL N >= ~400 at a per-SD effect >= ~0.3, or N ~ 5000 at effect ~ 0.12; GWAS N = 6000). Verified 2026-09-21 by simulation (r2 = 0.52, GWAS causal at SNP A + SNP B, eQTL causal at A, 20 replicates): coloc.abf did not blur to a middling PP.H3; it flipped per replicate between PP.H4 ~ 1 and PP.H3 ~ 1 (about 40% of replicates returned PP.H3 > 0.8 for a truly shared signal), depending on which of the two SNPs led in each dataset. Limited power does **not** trigger it: eQTL N <= 200, or small effects (0.05-0.12 at N <= 1000), gave PP.H1/PP.H0 dominance with PP.H3 <= 0.05. A well-powered single-causal locus at r2 0.5 is also decisive (PP.H4 = 1.0 shared, PP.H3 = 1.0 distinct), so r2 alone is not the trigger.

**Mechanism:** The single-causal-variant assumption forces each dataset onto one SNP. When the lead SNP differs between the GWAS (B) and the eQTL (A) because two nearby signals compete, the per-SNP Bayes factors of the two datasets peak at different SNPs and the posterior goes to H3 (distinct causal variants).

**Symptom:** LocusZoom looks convincing, but `result$summary['PP.H3.abf']` is high while PP.H4 is low, and the answer is unstable under resampling or a different eQTL panel. Because the outcome is near-binary (not intermediate), a decisive PP.H3 at moderate LD in a locus with 2+ GWAS signals is not proof of distinct causal variants: check with coloc.susie. Low PP.H4 with PP.H1 dominating is a power problem, not this failure.

**Fix:** Run coloc.susie (or eCAVIAR or PWCoCo) to allow multiple causal variants. If coloc.susie returns multiple credible sets with one pair showing PP.H4 > 0.75, this is real allelic heterogeneity not failure.

### Failure modes in reference files

- **coloc.susie -- LD reference mismatch** (`estimate_s_rss` lambda > 0.05, spurious credible sets): `references/coloc-susie.md`.
- **Lead-SNP swap and window bias** (PP.H4 changes when the window is re-centred): `references/lead-snp-swap.md`.

### coloc default p12 too liberal for trans-eQTL

**Trigger:** Applying `p12 = 1e-5` (the default) to a trans-eQTL / GWAS pair.

**Mechanism:** The default p12 was calibrated for cis-eQTL where biological proximity makes shared causality reasonable. For trans-eQTL, prior probability of shared causality is much lower; uniform p12 over-favours H4.

**Symptom:** PP.H4 > 0.8 reported, but sensitivity() reveals PP.H4 falls below 0.5 for p12 < 1e-5; replication in independent data fails.

**Fix:** Operational definition: "trans" = >5 Mb from TSS or different chromosome. Default p12=1e-5 over-favors H4 for trans (genome-rare biology). For trans: lower p12 to 5e-6 or 1e-6 AND raise PP.H4 threshold to >= 0.8 (compensate for higher FP risk). Cross-reference Vosa 2021 Nat Genet 53:1300 (eQTLGen trans) for empirical patterns.

### MHC / HLA + chr 8 inversion -- single-causal assumption breaks

**Trigger:** Locus within chr6:25-35 Mb (extended MHC, hg38), or chr8:8.1-11.9 Mb (chr 8 inversion, hg38).

**Mechanism:** The MHC contains classical HLA genes with extreme long-range LD (r2 > 0.5 over many Mb), multiple independent causal haplotypes, and structural variation. The chr 8p23.1 inversion similarly produces long-range LD across megabases of polymorphic inversion alleles. The single-causal-variant assumption is biologically wrong in both regions.

**Symptom:** coloc.abf almost always returns PP.H3 or fragmented PP across H1/H2/H3/H4 even when the underlying biology is well-established (e.g. HLA-DRB1 in autoimmune GWAS).

**Fix (MHC):** Use HLA-imputed classical alleles via SNP2HLA / HIBAG / HLA-TAPAS, then HLA-coloc (Butler-Laporte 2024 medRxiv) -- NOT coloc on SNPs in MHC. OR exclude MHC from genome-wide coloc and report HLA association at the haplotype/allele level. **Fix (chr 8 inversion):** Exclude chr8:8.1-11.9 Mb or pre-condition on inversion genotype before coloc. Never report a single coloc PP.H4 in either region without this caveat.

**Enforce programmatically, do not rely on remembering the coordinates:** call `flag_excluded_region()` (`scripts/flag_excluded_region.R`, already called by `scripts/coloc_abf.R` and `scripts/coloc_susie.R`) on the locus before running coloc.abf/coloc.susie, and abort/redirect rather than compute a naive PP.H4 when it returns non-NA.

### Underpowered eQTL (N < 200)

**Trigger:** Small eQTL discovery (e.g. tissue-specific bulk study, N < 200; per-cell-type sc-eQTL).

**Mechanism:** With low N, varbeta is large; the eQTL's per-SNP Bayes factors are flat; the joint likelihood concentrates on H0 or H1 (GWAS-only).

**Symptom:** PP.H0 or PP.H1 dominates; the eQTL panel shows visible signal but coloc cannot resolve causal vs noise.

**Fix:** Use eQTLGen (N ~ 31k whole-blood) or GTEx v8 (N ~ 70-700 per tissue) where possible. For rare cell types, accept the limitation and report the locus as underpowered rather than claim absence of colocalization.

| eQTL N | Coloc viability | Notes |
|--------|-----------------|-------|
| < 200 | Underpowered | PP.H1 dominant; flag |
| 200-500 | Cis only, modest | Single-tissue cis |
| 500-1000 | Good for cis | Most GTEx v8 tissues |
| >= 1000 | Well-powered | Trans accessible |
| >= 10000 | Meta (eQTLGen) | Cross-tissue / sc |

### Reference QTL panel choice

GTEx v8 (838 donors, 49 tissues, 2020) is the current PredictDB-supported standard. GTEx v10 (released 2024) has limited harmonisation and is not yet PredictDB-default. eQTLGen blood meta-eQTL (N ~ 31k) wins on sample size for blood cis-eQTL discovery, beating any single tissue on power. Always pin version in methods (e.g. "GTEx v8 MASHR-EUR, PredictDB release 2022-01").

## PP.H4 Threshold Framework

| Threshold | Use case | Source |
|-----------|----------|--------|
| 0.5 - 0.7 | Suggestive / pilot / hypothesis-generating | Giambartolomei 2014 original |
| >= 0.7 | Triangulation tier for TWAS / cis-MR / effector-gene cross-evidence | Open Targets Genetics common practice; cross-reference downstream skills |
| >= 0.75 | Open Targets Platform / eQTL Catalogue / FinnGen default screening threshold | Open Targets Genetics docs; Mountjoy 2021 Nat Genet 53:1527 |
| >= 0.80 | Most published colocalizations / standard publication tier | Wallace 2020 PLoS Genet 16:e1008720 |
| >= 0.90 | Stringent clinical / therapeutic-target prioritization | Reserved for high-confidence claims |
| >= 0.95 | Industry / regulatory drug-target submission grade | Internal pharma default |
| PP.H3 >= 0.80 | Confident distinct causal variants (negative coloc result) | Standard |
| PP.H4 / (PP.H3 + PP.H4) >= 0.9 | Conditional probability framing (some pipelines) | Foley 2021 |

**Operational rule:** Three operational tiers map onto the most common downstream uses: (a) **>= 0.7** when PP.H4 is one of several lines of triangulating evidence (TWAS + coloc, cis-MR + coloc, effector-gene multi-evidence) -- this is the threshold downstream skills (causal-genomics/transcriptome-wide-association, causal-genomics/mendelian-randomization cis-MR, causal-genomics/effector-gene-prioritization, causal-genomics/proteome-mr-drug-target) require; (b) **>= 0.8** for standard peer-reviewed publication as a stand-alone coloc claim (Wallace 2020); (c) **>= 0.95** for industry / clinical drug-target submission. Open Targets and FinnGen pipelines screen at >= 0.75 but downstream publication-grade coloc claims should clear >= 0.8 and triangulation claims >= 0.7. ALWAYS report PP.H3 alongside PP.H4 -- a locus with PP.H4 = 0.6, PP.H3 = 0.3 is qualitatively different from PP.H4 = 0.6, PP.H3 = 0.05 (the former is real ambiguity over single vs distinct causal; the latter is underpowered evidence). Run `coloc::sensitivity()` and report the p12 range over which PP.H4 stays above the threshold.

## Default Priors and the p12 Sensitivity Question

| Prior | Default | Interpretation | When to change |
|-------|---------|----------------|----------------|
| p1 | 1e-4 | Prob a random SNP is associated with trait 1 | Rarely changed |
| p2 | 1e-4 | Prob a random SNP is associated with trait 2 | Rarely changed |
| p12 | 1e-5 | Prob a random SNP is associated with both traits | Lower (5e-6 or 1e-6) for trans-eQTL or unrelated trait pairs; raise (5e-5) only with strong prior, e.g. molecular QTL in the same tissue as causal cell type |

The p12/p1 ratio (= 0.1 under defaults) is the prior odds of colocalization given a trait-1 association. Wallace 2020 (PLoS Genet 16:e1008720) showed default p12 = 1e-5 is too liberal for many real-world settings and recommended sensitivity analysis as standard practice. Pullin & Wallace 2025 (PLoS Genet 21:e1011697) extended this with variant-specific priors weighted by functional annotation.

### p12 Sensitivity Grid

| p12 grid point | Use case | Reporting rule |
|----------------|----------|-----------------|
| 1e-4 | Suggestive only / EUR cis-eQTL relaxed | PP.H4 here cannot support a publication claim |
| 1e-5 | Default for most cis-eQTL <-> GWAS pairs | Standard |
| 5e-6 | Conservative cis; default for trans-eQTL coloc | Recommended publication baseline |
| 1e-6 | Very conservative; trans coloc with weak prior | Required for cross-trait genome-rare coloc |

**Operational rule:** Require PP.H4 to remain above threshold across at least 3 adjacent grid points; report the lowest p12 at which PP.H4 >= 0.75. Use `coloc::sensitivity(result, rule = 'H4 > 0.75')` for the diagnostic plot.

Required reporting: PP.H4 at default priors + p12 range over which PP.H4 stays above threshold.

## Standard coloc.abf Pipeline

**Goal:** Test whether a single GWAS lead variant shares a causal variant with an eQTL gene's top signal at a defined locus.

**Approach:** Extract a 1 Mb window centred on the GWAS lead; harmonise alleles between datasets; format coloc input lists with `type` ('quant' or 'cc'), sample size `N`, and either `sdY` (quant) or `s` (cc); run coloc.abf; run sensitivity() over the p12 grid.

```bash
Rscript scripts/coloc_abf.R gwas.tsv eqtl.tsv --gwas-type cc --gwas-s 0.30 --gwas-n 50000     --eqtl-type quant --eqtl-sdy 1 --eqtl-n 500 --p12 5e-6 --out coloc_out
```

`scripts/coloc_abf.R` reads two harmonised tables (SNP, CHR, POS, BETA, SE; same SNP set and order; see `references/allele-harmonisation.md`), runs the MHC / chr 8 gate (`scripts/flag_excluded_region.R`, a hard `stop()`, run before any coloc call), builds the coloc input lists (`beta`, `varbeta`, `snp`, `position`, `type`, `N`, plus `s` for cc or `sdY` for quant), runs `coloc.abf` (p1 = p2 = 1e-4, conservative p12 = 5e-6) and `sensitivity(res, rule='H4 > 0.75')`, and writes `<out>_summary.tsv`, `<out>_snps.tsv` and `<out>_sensitivity.pdf`.

`sdY` semantics: when omitted, coloc estimates from `MAF` and `varbeta`; supplying `sdY=1` ASSUMES the trait is already standardised (eQTL with inverse-normal-transformed expression). Mismatch produces silently wrong Bayes factors -- the most common silent failure.

- For quantitative trait: leave `sdY=NULL` to estimate via `coloc:::sdY.est(varbeta, MAF, N)`. If CV of estimated sdY across SNPs > 0.5, varbeta/MAF are inconsistent -- coloc will silently miscalibrate Bayes factors.
- For inverse-normal-transformed expression: use `sdY = 1` (already standardized).
- Mismatch produces silently wrong PP -- most common silent failure.

`s` parameter for case-control (`type='cc'`):

- `s` = N_cases / N_total (NOT cases-per-control; NOT 0.5 default).
- For population-cohort case-control: `s` ~ disease prevalence in the cohort (~0.005 for rare disease).
- Wrong `s` does not error -- silently biases PP at extreme MAF.

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| `Error in coloc.abf: dataset must have N` | Forgot `N` in list, or `type` not set | Supply both; case-control also needs `s`; quant also needs `sdY` |
| PP.H3 dominant despite obvious visual overlap | 2+ causal in moderate LD breaks single-causal assumption | Run coloc.susie or eCAVIAR |
| `estimate_s_rss` lambda > 0.05 | LD reference does not match z-scores (ancestry / sample / build) | Use in-sample LD or ancestry-matched reference; do not proceed |
| PP.H4 unstable across p12 sensitivity grid | Borderline evidence; default priors not justified | Report the p12 range; lower priors for trans-eQTL; do not over-claim |
| coloc.susie returns NULL summary | No overlapping credible sets between traits | Genuine result (no shared signal) or both traits underpowered |
| Per-SNP betas have opposite signs but same magnitude across traits | Effect-allele mismatch | Run harmonisation; flip betas where A1/A2 swap; drop palindromic at high MAF |
| SMR significant + HEIDI p <= 0.05 | Linkage, not shared causal | Report as linkage; do not call colocalization |
| moloc all-share PPA collapses to ~0 | Different sample sizes / power across omics | Inspect per-omic effect sizes; consider HyPrColoc for cluster output |
| HyPrColoc trait cluster fragments | Underlying biology is multi-causal | Switch to coloc.susie at each suspected cluster centre |
| MHC PP.H4 close to 0 with strong visual signal | Long-range LD breaks single-causal | Use HLA-coloc or exclude MHC; never report standard coloc PP at MHC |

## Tool Install Notes

- **coloc**: CRAN. `install.packages('coloc')`. Bundles susieR dependency for >= 5.1.
- **susieR**: CRAN. `install.packages('susieR')`. >= 0.12.35 for `estimate_s_rss` and `kriging_rss`.
- **HyPrColoc**: GitHub only (never CRAN). `remotes::install_github('jrs95/hyprcoloc')`. Requires R >= 3.6. Its C++ (`align2.cpp`) fails to compile against RcppEigen 0.3.4.x (Eigen 3.4 `IndexedView` change: "cannot convert ... IndexedView ... to Scalar"); install `RcppEigen` 0.3.3.9.4 (Eigen 3.3.9) from the CRAN archive first. Built and run on Linux (R 4.4.1, g++ 15.2) that way; `knitr`/`rmarkdown` are listed in Imports but never used by the code.
- **SMR**: Pre-compiled binary from cnsgenomics.com/software/smr. Linux/Mac/Windows binaries; no R package.
- **eCAVIAR**: Compile from GitHub fhormoz/caviar; C++ source. CLI `eCAVIAR`. PAINTOR is the related multi-trait fine-mapping toolkit.
- **PWCoCo**: GitHub jwr-git/pwcoco. Compiled C++ CLI; can also be invoked from R via wrapper scripts.
- **SharePro_coloc**: GitHub only (no PyPI release). `git clone https://github.com/zhwm/SharePro_coloc` then `pip install -r requirements.txt`.
- **moloc**: GitHub clagiamba/moloc. R package; minimally updated since 2019, no CRAN release. R >= 3.5.
- **Plotting / data deps**: `install.packages(c('ggplot2', 'patchwork', 'data.table'))` -- ggplot2 + patchwork are used by `examples/regional_plots.R`.

## Reference Files

Read the file when the row or step it serves is in play; SKILL.md alone covers a coloc.abf run.

| File | Read when |
|------|-----------|
| `references/allele-harmonisation.md` | Before any coloc: merging GWAS and QTL stats, allele flips, strand, palindromic SNPs, `harmonise()` |
| `references/coloc-susie.md` | Two or more signals at the locus; building the signed LD matrix; the pipeline; LD-mismatch failure mode |
| `references/pwcoco.md` | Individual-level reference genotypes are available and GCTA-COJO finds >= 2 signals |
| `references/smr-heidi.md` | Pleiotropy vs linkage question; SMR vs coloc reconciliation; SMR command |
| `references/ecaviar-clpp.md` | Per-SNP CLPP under allelic heterogeneity; CLPP thresholds |
| `references/moloc.md` | Three to five omic layers at one locus |
| `references/hyprcoloc.md` | Many traits or tissues clustered at one locus |
| `references/lead-snp-swap.md` | The GWAS and QTL lead SNPs differ; diagnosing window-centring bias |
| `references/reporting.md` | Writing up a result: reviewer pushback responses and the reporting template |

## Scripts

| Script | Use |
|--------|-----|
| `scripts/harmonise.R` | Harmonise two summary-stat tables (CLI or `source()` for `harmonise()`) |
| `scripts/flag_excluded_region.R` | MHC / chr 8 inversion gate; sourced by the two coloc scripts |
| `scripts/coloc_abf.R` | Gate, `coloc.abf`, p12 sensitivity for one locus |
| `scripts/coloc_susie.R` | Gate, LD-consistency check, `runsusie` x2, `coloc.susie` |

Runnable demos on simulated data are in `examples/`.

## References

- Giambartolomei C et al 2014 PLoS Genet 10:e1004383 (coloc.abf)
- Wallace C 2020 PLoS Genet 16:e1008720 (prior elicitation; relaxing the single-causal-variant assumption; default-prior sensitivity)
- Pullin JM & Wallace C 2025 PLoS Genet 21:e1011697 (variant-specific priors)
- Wallace C 2021 PLoS Genet 17:e1009440 (coloc.susie; multiple causal variants)
- Zhu Z et al 2016 Nat Genet 48:481 (SMR + HEIDI)
- Hormozdiari F et al 2016 AJHG 99:1245 (eCAVIAR / CLPP)
- Giambartolomei C et al 2018 Bioinformatics 34:2538 (moloc)
- Foley CN et al 2021 Nat Commun 12:764 (HyPrColoc)
- Robinson JW et al 2022 bioRxiv 2022.08.08.503158 (PWCoCo)
- Zhang W et al 2024 Bioinformatics 40:btae295 (SharePro_coloc)
- Butler-Laporte G et al 2024 medRxiv 2024.11.05.24316783 (hlacoloc)
- Mountjoy E et al 2021 Nat Genet 53:1527 (Open Targets Genetics colocalization pipeline)
- Vosa U et al 2021 Nat Genet 53:1300 (eQTLGen, N ~ 31,684 whole blood)
- GTEx Consortium 2020 Science 369:1318 (GTEx v8 multi-tissue eQTL)

## Related Skills

- causal-genomics/mendelian-randomization - Causal effect estimation from coloc-validated SNPs
- causal-genomics/fine-mapping - SuSiE / FINEMAP / CAVIAR credible sets feeding coloc.susie; LD construction protocol cross-ref
- causal-genomics/mediation-analysis - Downstream causal mediation given coloc shared causal variants
- causal-genomics/pleiotropy-detection - Distinguishing horizontal pleiotropy from shared causality
- causal-genomics/transcriptome-wide-association - TWAS / PrediXcan / FOCUS gene-level prioritization complementary to coloc
- causal-genomics/proteome-mr-drug-target - pQTL coloc + MR for drug-target prioritization
- causal-genomics/effector-gene-prioritization - Locus-to-gene with coloc, ABC, V2G integration
- population-genetics/association-testing - GWAS summary statistic generation and locus extraction
- population-genetics/linkage-disequilibrium - LD reference panels for coloc.susie and PWCoCo
- variant-calling/variant-annotation - Functional annotation for variant-specific priors
- variant-calling/filtering-best-practices - Pre-coloc QC for summary stats
- differential-expression/deseq2-basics - Generating eQTL / molecular QTL counts
- single-cell/scatac-analysis - Per-cell-type chromatin context for coloc interpretation
- workflows/gwas-pipeline - Upstream GWAS analysis producing coloc input
- data-visualization/ggplot2-fundamentals - Regional and LocusCompare plot construction
