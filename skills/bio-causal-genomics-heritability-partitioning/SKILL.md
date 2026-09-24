---
name: bio-causal-genomics-heritability-partitioning
description: Estimates SNP heritability and partitions it across functional annotations, cell types, and loci from GWAS summary statistics or individual-level genotypes. Implements LDSC, stratified LDSC with the baseline-LD model, Finucane 2018 cell-type prioritization, LDAK SumHer, HDL, HESS local heritability, BOLT-REML, GCTA-GREML, graphREML, and Popcorn cross-population genetic correlation. Use when computing total h2_SNP from summary stats, partitioning heritability across functional categories, prioritizing trait-relevant tissues or cell types from ENCODE/Roadmap chromatin marks, reconciling LDSC vs LDAK enrichment estimates, computing local heritability with HESS, estimating genetic correlation between traits, or producing publication-grade enrichment with calibrated sensitivity to model assumptions.
tool_type: mixed
primary_tool: ldsc
license: MIT
author: GPTomics
---

## Version Compatibility

Reference examples tested with: LDSC v3.0.1 (`CBIIT/ldsc`, `main` branch, commit `1f09cf0c`, Python
3.9+ -- see "Tool Install Notes"), LDAK 6.3, BOLT-LMM 2.5 (see caveat below), GCTA 1.94.1 (verified
2026-09-21, see "Per-Method Reference Files"), HESS 0.5.4-beta (verified 2026-09-21 with the Python 3
patches in `references/hess-local-h2.md` -- upstream has no separate "Python 3 branch"), HDL 1.4.0+ (R;
GitHub `zhenin/HDL`; blocked here on a ~5 GB reference-panel download, not attempted), Popcorn 1.1
(verified 2026-09-21, see `references/popcorn-transancestry.md`), baselineLD_v2.2 annotations
(alkesgroup.broadinstitute.org/LDSCORE).

**BOLT-LMM download URL was stale and is now fixed**: `alkesgroup.broadinstitute.org/BOLT-LMM/downloads/BOLT-LMM_v2.4.1.tar.gz`
404s -- the current release moved to `.../downloads/BOLT-LMM_v2.5.tar.gz`, with older versions under
`.../downloads/old/`. Checked 2026-09-21. **Also found (2026-09-21): the v2.5 Linux binary segfaults
immediately on this machine's WSL environment, even on `bolt --help`** -- not confirmed whether this is
a CPU-feature or glibc incompatibility specific to this host; building from the bundled `src/` was not
attempted (needs a C++ toolchain and Intel MKL). Treat BOLT-REML as unverified on this environment;
verify the binary with `bolt --help` before relying on it elsewhere.

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `python -c 'import <module>; help(<module>)'`
- R: `packageVersion('<pkg>')` then `?function_name`
- CLI: `<tool> --version` then `<tool> --help`

**LDSC fork history (checked 2026-09-17):** `bulik/ldsc`'s own README was updated 2026-01-16 to say
it is superseded by `CBIIT/ldsc` (NIH Center for Biomedical Informatics and Information Technology),
a maintained Python 3 port with a 2026-08-10 commit and an accompanying preprint. Earlier advice to
use `abdenlab/ldsc-python3` or `belowlab/ldsc` is retracted: both were verified broken on `--h2`,
`--rg`, and/or `--h2-cts` in the exact pinned version (wrong internal dispatch names, a numpy>=2
crash, an `n_chr`/`num` keyword mismatch). Use `CBIIT/ldsc` per "Tool Install Notes" below; `--h2` and
`--rg` run unmodified, `--h2-cts` needs one documented one-line patch (same section). If code throws
ImportError, AttributeError, or a "category not found" error in the LD score file, introspect the
installed binary and the actual LD-score column headers rather than retrying.

## Scope

h2 and PRS outputs are population-level research statistics (variance explained across a cohort), not
an individual diagnostic or prescriptive tool -- never use them to tell a specific person their odds of
developing a disease or what treatment to start; redirect that request to a clinician or genetic
counselor. If asked for a heritability or enrichment number with no data supplied, do not invent one;
either run the real pipeline on data the user provides, or offer a clearly-cited published estimate
from the literature and label it as external, not computed.

# Heritability Partitioning

**"Estimate SNP heritability and partition it across functional categories, cell types, and loci"** -> Decompose `h2_SNP` from GWAS summary statistics (or individual-level genotypes) into contributions from baseline annotations (coding, conserved, regulatory), tissue-specific chromatin marks, and per-locus components, then reconcile model-dependent enrichment estimates across LDSC and LDAK. Tool choice is a decision about the **regime** (summary-stat vs individual-level; one-trait vs two-trait genetic correlation; total vs partitioned vs local) and the **model assumption** about how per-SNP heritability scales with LD, MAF, and functional annotation (GCTA model vs LDAK-Thin vs baseline-LD).

- CLI (h2 from sumstats, EUR): `ldsc.py --h2 trait.sumstats.gz --ref-ld-chr eur_w_ld_chr/ --w-ld-chr eur_w_ld_chr/ --out h2`
- CLI (functional partitioning): `ldsc.py --h2 trait.sumstats.gz --ref-ld-chr baselineLD.,<annot>. --frqfile-chr 1000G.EUR.QC. --w-ld-chr weights. --overlap-annot --print-coefficients --out part`
- CLI (cell-type prioritization, Finucane 2018): `ldsc.py --h2-cts trait.sumstats.gz --ref-ld-chr-cts <cts_file>.ldcts --w-ld-chr weights. --out cts`
- CLI (cross-trait rg): `ldsc.py --rg t1.sumstats.gz,t2.sumstats.gz --ref-ld-chr eur_w_ld_chr/ --w-ld-chr eur_w_ld_chr/ --out rg`
- CLI (LDAK alternative): `ldak --sum-hers <out> --summary trait.txt --tagfile ldak.thin.<build>.tagging --check-sums NO`
- R (HDL): `HDL::HDL.rg(gwas1.df, gwas2.df, LD.path = 'UKB_array_SVD_eigen90_extraction')`
- CLI (local h2): HESS step1 `hess.py --local-hsqg trait.sumstats.gz --chrom <chr> --bfile <ref> --partition <part>.bed --out hess_<chr>`

## Statistical Model Taxonomy

| Method | Heritability model | Input | Output | Strength | Fails when |
|--------|--------------------|-------|--------|----------|------------|
| LDSC h2 (Bulik-Sullivan 2015 Nat Genet 47:291) | GCTA model: per-SNP h2 proportional to LD score | Sumstats + ancestry-matched LD scores | h2 estimate + intercept + ratio | Standard for sumstats; fast; calibrated EUR LD scores | N too low (mean chi-square < 1.02); non-EUR ancestry with EUR LD scores; population stratification not captured |
| Stratified LDSC / S-LDSC (Finucane 2015 Nat Genet 47:1228) | GCTA model with per-annotation tau coefficients | Sumstats + baseline + custom annotations | Per-annotation enrichment + tau | Reference functional partitioning method | Highly collinear annotations inflate per-tau SE; small annotation (<0.5% genome) underpowered |
| Baseline-LD model (Gazal 2017 Nat Genet 49:1421) | Adds LD-related and MAF-dependent annotations to baseline | Sumstats + `baselineLD_v2.2.` | Enrichment robust to LD-MAF confounding | Modern S-LDSC default; calibrates LD/MAF dependence | EUR-only baseline-LD v2.2 must NOT be used on non-EUR GWAS; use baseline-LD-X / S-LDXR (Shi H & Gazal S et al 2021 Nat Commun 12:1098) instead |
| Cell-type S-LDSC (Finucane 2018 Nat Genet 50:621) | Per-cell-type annotation marginal to baseline | Sumstats + cell-type chromatin annotations (.ldcts) | Per-cell-type p-value | Tissue / cell-type prioritization; published per-tissue ldcts files | Sample size small (mean chi-square < 1.02); annotation overlaps strongly with baseline |
| HDL (Ning 2020 Nat Genet 52:859) | Genome-wide eigen-decomposition likelihood | Sumstats + HDL reference panel (UKB N=336k) | h2 and rg with ~60% lower variance than LDSC | Equivalent to ~2.5x sample size for h2 / rg | Sample overlap > 5% biases the likelihood; only EUR HDL reference panel publicly available; no non-EUR HDL eigen-reference exists as of 2026 -- for non-EUR fall back to ancestry-matched cross-trait LDSC (intercept absorbs overlap; rg unbiased) |
| LDAK SumHer (Speed 2019 Nat Genet 51:277) | LDAK-Thin model: per-SNP h2 reweighted by MAF + LD | Sumstats + LDAK-Thin tagging file | h2 + enrichment | Alternative to LDSC; often better-fitting per cross-validation per Speed 2017 Nat Genet 49:986 | Tagging file must match build / ancestry; non-EUR support limited |
| HESS (Shi 2016 AJHG 99:139; Shi 2017 AJHG 101:737) | Per-locus h2 via quadratic form on LD-projected effect estimates | Sumstats + LD reference + locus partition (LDetect) | Per-locus h2 + bivariate local rg | Locus-level resolution; identifies high-h2 loci for follow-up | Locus has < 1000 SNPs; LD reference must be in-sample or matched |
| BOLT-REML (Loh 2015 Nat Genet 47:1385) | Bayesian REML; multi-component variance | Individual-level genotypes (PLINK BED) + phenotype | h2 + per-component partition | Biobank-scale (N=500k feasible); more precise than LDSC at high N | Needs individual-level data; assumes Gaussian residual; not for case-control < 5% prevalence without LMM-BOLT |
| GCTA-GREML (Yang 2011 AJHG 88:76) | GRM-based REML on individual genotypes | GRM (PLINK BED) + phenotype + covariates | h2 + SE | Gold standard for individual-level data; PCGC for case-control | N <= 5000 has wide SE; case-control needs PCGC-S correction; population stratification leaks into h2 |
| graphREML (Li H et al 2024 medRxiv 2024.11.04.24316716; published Nat Genet 2026) | Sumstat REML using LDGM graph | Sumstats + LDGM reference | h2 + functional partition | Use when an LDGM reference exists for the ancestry AND runtime matters at biobank N (> 200k); pre-built LDGM panels currently cover EUR + EAS | Newer (2024); reference panel availability evolving; non-EUR/EAS ancestries lack pre-built LDGM |
| Popcorn (Brown 2016 AJHG 99:76) | Trans-ancestry genetic correlation under MAF-LD model | Sumstats + cross-population LD scores | rg trans-ancestry + h2 per population | Distinguishes shared-causal vs ancestry-specific signal | Effective N per population must be > 5000; small non-EUR cohorts unstable |
| Cross-model reconciliation (Gazal 2019 Nat Genet 51:1202) | Joint LDSC + LDAK comparison framework | Both LDSC and LDAK outputs | Enrichment model-comparison (Gazal: baseline-LD better-calibrated; LDAK developers dispute) | Quantifies model-dependent component of enrichment | Methodological / reporting practice, not a separate primary estimator |

Methodology evolves; benchmark consensus shifts. Verify against current Yengo 2022 *Nat Methods*, Speed D, Holmes J & Balding DJ 2020 Nat Genet 52:458 (heritability model comparison), and the alkesgroup/Price Lab tutorial (LDSC) before locking method as primary. Per-SNP-h2 model choice is an open debate; report both LDSC and LDAK SumHer when a claim depends on model assumption.

## Decision Tree by Scenario

| Scenario | Recommended | Why |
|----------|-------------|-----|
| Total h2 from sumstats, EUR GWAS, N > 50k | LDSC h2 | Standard, fast, well-calibrated against EUR reference |
| Partitioned h2 by functional category | S-LDSC with baseline-LD_v2.2 | Default functional partitioning; LD/MAF-robust |
| Tissue / cell-type prioritization | S-LDSC --h2-cts with ENCODE/Roadmap or scATAC ldcts | Designed for this; per-tissue Bonferroni-controlled; workflow, `.ldcts` format in `references/cell-type-prioritization.md` |
| Two-trait genetic correlation from sumstats, no overlap | HDL (primary) + cross-trait LDSC (secondary) | HDL ~60% lower variance; LDSC robust under any overlap; HDL in `references/hdl-genetic-correlation.md`, LDSC rg in `references/cross-trait-ldsc.md` |
| Two-trait rg with sample overlap > 5% | Cross-trait LDSC | HDL biased by overlap; LDSC intercept absorbs overlap; how to run and read `gcov_int` in `references/cross-trait-ldsc.md` |
| Individual-level biobank h2, N > 100k | BOLT-REML | Better precision; multi-component partition |
| Smaller individual-level cohort, N 5-50k | GCTA-GREML | Gold-standard REML; PCGC if case-control < 20% prevalence; pipeline in `references/gcta-greml.md` |
| Local heritability and bivariate local rg | HESS | Per-locus resolution; identifies hotspots for follow-up; pipeline in `references/hess-local-h2.md` |
| Trans-ancestry rg / cross-population h2 | Popcorn | Designed for trans-ethnic; LD scores per population; pipeline in `references/popcorn-transancestry.md` |
| Functional enrichment claim depends on model | Report BOTH LDSC and LDAK SumHer | Per Gazal 2019; model-dependence is real; LDAK pipeline in `references/ldak-sumher.md` |
| Case-control GWAS with low prevalence | LDSC on liability scale (--samp-prev --pop-prev) | Observed-scale h2 understates liability-scale truth |
| Single-cell ATAC cell-type prioritization | S-LDSC with per-cluster ATAC peaks as annotations | Cross-reference atac-seq/single-cell-atac for peak generation |

## Per-Method Reference Files

Load only the file for the method in use (each holds the pipeline, install line, failure modes and errors for that method):

| Method | File |
|--------|------|
| LDAK SumHer (`--calc-tagging`, `--sum-hers`) | `references/ldak-sumher.md` |
| HESS local h2 | `references/hess-local-h2.md` |
| HDL genetic correlation (R) | `references/hdl-genetic-correlation.md` |
| Cell-type / tissue prioritization (`--h2-cts`, `.ldcts` manifest, Bonferroni) | `references/cell-type-prioritization.md` |
| Cross-trait LDSC rg (`--rg`, `gcov_int` under sample overlap) | `references/cross-trait-ldsc.md` |
| GCTA-GREML (individual-level h2) | `references/gcta-greml.md` |
| Popcorn (trans-ancestry rg) | `references/popcorn-transancestry.md` |
| Failure modes: intercept misread, non-EUR LD scores, collinear annotations | `references/ldsc-failure-modes.md` |
| Runtime and hardware per method | `references/computational-footprint.md` |

Total h2, partitioned S-LDSC and the case-control liability rule stay in this file. BOLT-REML and
graphREML carry no pipeline here beyond the taxonomy row, decision tree and install notes (BOLT-REML's
binary is unverified on this Skill's own environment -- see "Version Compatibility"; graphREML needs an
LDGM reference this environment does not have).

## LDSC Intercept Interpretation (Postdoc Nuance)

The LDSC intercept is widely misinterpreted as a "confounding score". The correct interpretation:

- Intercept = 1 indicates no inflation from population structure, cryptic relatedness, or sample overlap (idealised)
- Intercept > 1 indicates SOME source of inflation, BUT polygenic background can elevate the intercept too: at very high N, polygenic signal can lift the intercept modestly without stratification
- The **ratio** statistic `ratio = (intercept - 1) / (mean_chi2 - 1)` is the fraction of inflation attributable to non-polygenic sources; a ratio of 0 means all inflation is polygenic, a ratio near 1 means most of it is stratification or overlap
- Bulik-Sullivan 2015 recommends interpreting intercept jointly with mean chi-square; intercept of 1.05 is innocuous if mean chi-square is 1.5 (ratio = 0.1) but worrying if mean chi-square is 1.05 (ratio = 1.0)

**Operational rule:** Always report intercept, mean chi-square, and ratio together. Do not interpret intercept in isolation. For sample-overlap diagnosis between two GWAS, use bivariate LDSC intercept, not univariate.

| Intercept | Ratio | Interpretation |
|-----------|-------|----------------|
| ~1.0 | ~0 | No inflation; h2 trustworthy |
| 1.0 - 1.1 | < 0.2 | Mostly polygenic; h2 trustworthy |
| 1.1 - 1.3 | 0.2 - 0.5 | Mild inflation; investigate population structure / sample overlap |
| 1.3 - 1.5 | 0.5 - 0.8 | Substantial inflation; report ratio jointly; consider re-genotype-QC |
| > 1.5 | -- | Re-run after PC adjustment or genomic control; do not interpret h2 |

**Intercept > 1.5 troubleshooting ladder** (work in order; stop when source is found): (a) per-cohort PC adjustment was insufficient; refit GWAS with more PCs (10-20) or per-cohort separately, (b) cryptic relatedness in the GWAS cohort -- run `king --related` and remove pairs with kinship > 0.05 (or 0.0884 for second-degree), (c) case-control matching imbalance -- check case/control PCs separately, (d) sample-overlap with one of the contributing cohorts (especially in meta-analysis) -- check bivariate intercepts pairwise, (e) if biobank-internal, recompute the GWAS with sample-level relatedness exclusion before LDSC.

## LDSC vs LDAK Reconciliation

LDSC (GCTA model) and LDAK SumHer (LDAK-Thin model) make different assumptions about how per-SNP heritability scales with LD and MAF:

- **GCTA model (LDSC default):** per-SNP h2 inversely proportional to local LD score; high-LD SNPs tag many causal variants
- **LDAK-Thin (Speed 2017 Nat Genet 49:986-992 introduced the LDAK model; the Thin model was formalized in Speed 2020 Nat Genet 52:458):** per-SNP h2 weighted by MAF and inversely by LD with empirical exponents; less weight on common high-LD SNPs

These give systematically different functional enrichment estimates. Speed 2019 reported that LDAK-Thin often gives **lower** conserved-region enrichment than baseline LDSC; conversely LDSC can over-attribute h2 to coding/conserved regions because the GCTA prior couples LD to causality. Gazal 2019 Nat Genet 51:1202-1204 argues the baseline-LD S-LDSC model is better-calibrated (higher model likelihood) and that LDAK/SumHer's lower functional-enrichment estimates are downward-biased; the LDAK developers dispute this (Speed 2020 Nat Genet 52:458), so report both models and flag the model-dependence.

**Operational rule:** Whenever functional enrichment is the primary claim (e.g. "h2 is enriched in tissue T by N-fold"), report enrichment from BOTH LDSC and LDAK. Flag the model assumption. If LDSC and LDAK disagree by > 2x (e.g. LDSC conserved-region enrichment 25x vs LDAK 10x; both are model-internally consistent and the data alone cannot pick one), treat the claim as model-dependent and cite Gazal 2019. Report LDSC primary + LDAK confirmatory, emphasize directional agreement over magnitude, and never pick the model that gives the desired answer. For non-enrichment claims (total h2, rg between two traits), the model dependence is smaller and LDSC alone is acceptable.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| mean chi-square > 1.02 | LDSC wiki / Bulik-Sullivan 2015 | Below this, LDSC h2 estimate has huge SE; need N proportional to 1 / h2 |
| Total h2 must lie in [0, 1] on both observed and liability scale | Definition of h2_SNP | A point estimate outside [0, 1] is never a face-value result. If the 95% CI (estimate +/- 1.96 SE) still excludes [0, 1] (e.g. liability h2 1.43 (0.10)), it is an artifact: check `--samp-prev`/`--pop-prev` assignment, munged N and the LD reference, then report it as a data or model failure. If the CI overlaps the bound, report it as noise around the boundary |
| h2 SE < 0.02 | LDSC convention | Below this SE, h2 estimate is interpretable; above, treat as exploratory |
| h2 SE >= 0.02 OR mean chi-square < 1.02 | LDSC convention | Estimate unreliable; N >= 50k is typical noise floor for h2 ~ 0.1 (scales as ~1/h2) |
| LDSC intercept / ratio bins | Bulik-Sullivan 2015 | See the intercept/ratio table under "LDSC Intercept Interpretation" |
| Stratified LDSC enrichment p < 0.05 / N_annot | Finucane 2015 | Bonferroni across annotations in the baseline-LD model |
| S-LDSC cell-type p < 2.5e-4 | Finucane 2018; ~200 tissues | Bonferroni for tissue prioritization |
| HESS h2 per locus needs >= 1000 SNPs | Shi 2017 AJHG 101:737 | Quadratic form unstable below this density |
| HDL sample overlap < 5% | Ning 2020 Nat Genet 52:859 | Above 5%, HDL likelihood is biased |
| LDAK tagging file build match | Speed 2019 | hg19 vs hg38 tagging files non-interchangeable |
| Annotation > 0.5% of genome | Finucane 2015 | Smaller categories underpowered for tau estimation |
| Effective N > 5000 per population (Popcorn) | Brown 2016 AJHG 99:76 | Below this, trans-ancestry rg has very wide CI |

## LDSC Standard Workflow

**Goal:** Compute total h2 plus partitioned heritability across functional categories from EUR GWAS summary statistics.

**Approach:** Munge sumstats to LDSC format -> run --h2 for total -> run --h2 with --ref-ld-chr including baseline annotations -> --overlap-annot for enrichment p-values -> --print-coefficients for per-annotation tau.

Run it as `bash examples/ldsc_partitioned_h2.sh <gwas.tsv> <trait_prefix> [samp_prev] [pop_prev]`: step 1 `munge_sumstats.py` to `<trait>.sumstats.gz` (SNP, A1, A2, Z, N), step 2 total h2 (intercept, ratio, mean chi2; liability scale when both prevalences are given), step 3 partitioned h2 against baseline-LD v2.2 (`--overlap-annot --print-coefficients`). Steps 4-5 of the same script are the cell-type and rg runs. It expects the reference bundles from "Tool Install Notes" under `$LDSC_REF`.

Ancestry-matched LD scores (EAS, AFR, AMR) are available at alkesgroup.broadinstitute.org/LDSCORE; do NOT apply EUR LD scores to non-EUR GWAS.

## Per-Method Failure Modes

LDSC intercept misinterpretation, non-EUR ancestry with EUR LD scores and collinear annotations: read `references/ldsc-failure-modes.md`.

### Case-control LDSC observed vs liability scale

**Trigger:** Reporting LDSC h2 from a case-control GWAS on the observed (0/1) scale.

**Mechanism:** Observed-scale h2 depends on case fraction in the GWAS sample, not population prevalence; comparisons across studies require liability-scale transformation.

**Symptom:** h2 looks tiny (0.02) for a known-heritable disease; differs across studies with different case fractions.

**Fix:** Always supply `--samp-prev <case_fraction>` and `--pop-prev <population_lifetime_prevalence>` to LDSC; report h2 on liability scale. Without these flags, LDSC defaults to observed scale. Lee 2011 (AJHG 88:294) conversion: `h2_liab = h2_obs * (K(1-K))^2 / (P(1-P) * z^2)` (numerator is K^2(1-K)^2), where K = population prevalence, P = sample case proportion, z = standard-normal density at the quantile (1-K). LDSC applies this via `--samp-prev/--pop-prev`; verify K and P are assigned correctly (K is population, P is sample). Skipping the conversion typically yields a 2-10x underestimate vs the liability-scale truth.

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| LDSC h2 negative | Trait truly null heritable; or LD scores ancestry-mismatched | Verify mean chi-square > 1.02; switch to ancestry-matched LD; report as null |
| `munge_sumstats.py` drops > 50% of SNPs | Allele mismatch with HapMap3 SNPlist; or A1/A2 swapped | Inspect drop reasons; pre-harmonize alleles; check rsID format |
| LDSC intercept > 1.5 | Population stratification, cryptic relatedness, or extensive sample overlap | Report jointly with ratio; investigate per-cohort PC adjustment |
| Stratified LDSC enrichment > 100x | Tiny annotation (<0.1% genome) underpowered | Set annotation size floor >= 0.5%; report joint enrichment of larger composite category |
| `--samp-prev/--pop-prev` not supplied for case-control | LDSC defaults to observed scale | Always supply both for case-control; report liability-scale h2 |
| Cell-type S-LDSC: no tissue p < 2.5e-4 | Trait underpowered for tissue prioritization (mean chi-square low) | Pool with related traits via MTAG; or report no tissue distinguishable |
| BOLT-REML "matrix not positive definite" | GRM has duplicated individuals or extreme relatedness | Pre-filter to unrelated < 0.05 kinship; or use REML method-of-moments instead |
| `munge_sumstats.py`: missing N column | Per-SNP N column not supplied | Pass `--N <Ntot>` (numeric total) OR `--N-col N` (column name); pick one |
| `munge_sumstats.py`: all SNPs drop silently | A1/A2 flipped relative to HM3 reference (`--merge-alleles`) | Verify allele coding matches `w_hm3.snplist`; pre-harmonize or swap A1 <-> A2 |
| `munge_sumstats.py` drops multi-allelic SNPs | Expected behavior (LDSC requires biallelic) | Document the drop count in methods; not a fix |
| Sign of Z reversed across studies | `--signed-sumstats BETA,0` vs `--signed-sumstats Z,0` confusion | Use `BETA,0` for additive effect sign convention, `Z,0` for Z-score; never both; verify direction with a known sentinel SNP |

## Tool Install Notes

```bash
# LDSC: use CBIIT/ldsc (NIH CBIIT's maintained Python 3 fork; bulik/ldsc's own README
# points here as of 2026-01-16). Checked on commit 1f09cf0 (2026-08-10), Python 3.9.
git clone https://github.com/CBIIT/ldsc.git
cd ldsc
conda create --name ldsc python=3.9 -y && conda activate ldsc   # or: micromamba create -n ldsc -c conda-forge -c bioconda python=3.9 bitarray=2 pybedtools=0.10.0 -y
pip install numpy==1.21.5 pandas==1.3.3 scipy==1.7.3   # environment3.yml's pins; conda-installed numpy/pandas float otherwise
./ldsc.py -h   # verify: prints the full flag list, no traceback

# --h2 and --rg run as-is against this clone. --h2-cts needs one line patched first:
#   ldscore/sumstats.py ~line 285: ref_ld_cts = ... .loc[:,1:])  ->  ... .iloc[:,1:])
# (a leftover from replacing pandas' removed .ix[:,1:] with the label-based .loc instead of
# the positional .iloc; pandas >=1.0 raises TypeError on the unpatched line.) Without the
# patch, --h2-cts fails with "cannot do slice indexing on Index with these indexers [1] of
# type int". LD score files must be gzipped (`.l2.ldscore.gz`) -- ldscore() hardcodes that
# suffix; every official reference bundle already ships gzipped.

# Pre-computed reference resources (one-time download)
# EUR LD scores (HapMap3 SNPs)
wget https://alkesgroup.broadinstitute.org/LDSCORE/eur_w_ld_chr.tar.bz2
# Baseline-LD v2.2 (Gazal 2017)
wget https://alkesgroup.broadinstitute.org/LDSCORE/1000G_Phase3_baselineLD_v2.2_ldscores.tgz
# 1000G EUR frequency files
wget https://alkesgroup.broadinstitute.org/LDSCORE/1000G_Phase3_frq.tgz
# Weights
wget https://alkesgroup.broadinstitute.org/LDSCORE/1000G_Phase3_weights_hm3_no_MHC.tgz
# Multi-tissue chromatin ldcts (Finucane 2018)
wget https://alkesgroup.broadinstitute.org/LDSCORE/Multi_tissue_chromatin_1000Gv3_ldscores.tgz
```

Before pointing LDSC at real GWAS data, confirm the install works: `LDSC_DIR=./ldsc bash examples/smoke_test_ldsc.sh`
runs `--h2`, `--rg`, and `--h2-cts` against the simulated fixtures in the clone's own `test/` directory
and prints real, checkable regression output for each.

LDAK, HESS and HDL install notes live in their `references/` files (see "Per-Method Reference Files").

```bash
# BOLT-LMM / BOLT-REML (checked 2026-09-21: the old v2.4.1 URL below 404s; current release is v2.5,
# older versions moved under .../downloads/old/)
wget https://storage.googleapis.com/broad-alkesgroup-public/BOLT-LMM/downloads/BOLT-LMM_v2.5.tar.gz
# Segfaults immediately (even --help) on this Skill's own verification environment (WSL, 2026-09-21) --
# not yet root-caused; verify with `bolt --help` before trusting the binary on a new machine.
```

GCTA-GREML and Popcorn each have their own pipeline, install and worked example: see "Per-Method
Reference Files" below (`references/gcta-greml.md`, `references/popcorn-transancestry.md`).

## References

- Bulik-Sullivan BK et al 2015 Nat Genet 47:291 (LDSC)
- Bulik-Sullivan BK et al 2015 Nat Genet 47:1236 (cross-trait LDSC genetic correlation)
- Finucane HK et al 2015 Nat Genet 47:1228 (stratified LDSC baseline)
- Gazal S et al 2017 Nat Genet 49:1421 (baseline-LD model)
- Shi H & Gazal S et al 2021 Nat Commun 12:1098 (S-LDXR / baseline-LD-X cross-population)
- Finucane HK et al 2018 Nat Genet 50:621 (cell-type S-LDSC)
- Speed D et al 2017 Nat Genet 49:986-992 (LDAK model; cross-validation model comparison)
- Speed D, Holmes J & Balding DJ 2020 Nat Genet 52:458 (LDAK-Thin model; heritability model comparison)
- Speed D et al 2019 Nat Genet 51:277 (SumHer)
- Gazal S et al 2019 Nat Genet 51:1202-1204 (LDSC vs LDAK reconciliation, model-dependence)
- Hou K et al 2019 Nat Genet 51:1244 (heritability accuracy and h2 model comparison)
- Ning Z et al 2020 Nat Genet 52:859 (HDL)
- Shi H et al 2016 AJHG 99:139 (HESS univariate)
- Shi H et al 2017 AJHG 101:737 (HESS bivariate local rg)
- Loh PR et al 2015 Nat Genet 47:1385 (BOLT-REML)
- Yang J et al 2011 AJHG 88:76 (GCTA-GREML)
- Brown BC et al 2016 AJHG 99:76 (Popcorn trans-ancestry rg)
- Berisa T & Pickrell JK 2016 Bioinformatics 32:283 (LDetect locus partition)
- Li H, Kamath T, Mazumder R, Lin X, O'Connor LJ 2024 medRxiv 2024.11.04.24316716 (graphREML; published Nat Genet 2026)

Related Skills are listed in `usage-guide.md`.
