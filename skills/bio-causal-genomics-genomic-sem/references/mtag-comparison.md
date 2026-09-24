# MTAG comparison and reconciliation (reference for genomic-sem SKILL.md)

Running MTAG beside GenomicSEM, and what to do when the two disagree. Read it when `SKILL.md` "Reference Files" points here. The section names below refer to `SKILL.md`.

## MTAG Comparison

**Goal:** Cross-check GenomicSEM common-factor results against MTAG per-trait shrunk z-scores.

**Approach:** Run MTAG CLI on the same input sumstats; compare top hits with GenomicSEM factor hits. Report MaxFDR.

Run it with `examples/mtag_pipeline.sh`: the same `mtag.py` call (`--n_min 0`, signed Z; `--use_beta_se` is disabled upstream and raises a RuntimeError since Dec 2021), the `max ?fdr` grep of the MTAG log (matches both the section header and the `Max FDR of Trait` value lines), and a genome-wide-significant hit count per trait. Each per-trait MTAG file is `<out>_trait_<k>.txt`.

If MaxFDR > 0.05 for any trait, MTAG results for that trait are unreliable; GenomicSEM with Q_SNP filtering is the more defensible report.

## Reconciliation: When GenomicSEM and MTAG Disagree

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| GenomicSEM factor SNP sig, MTAG sig for all traits | Genuine common-factor SNP | Report; high confidence |
| GenomicSEM factor SNP sig, MTAG sig in only 1 trait | Q_SNP heterogeneity likely; one-trait-dominant | Check Q_SNP; if sig, this is NOT a factor SNP |
| MTAG sig, GenomicSEM factor null, Q_SNP sig | Trait-specific SNP captured by MTAG shrinkage | Report as trait-specific, not common-factor |
| Both null but per-trait univariate sig | Power loss from multivariate parameterization | Re-check sample overlap V matrix |
| GenomicSEM and MTAG both sig but opposite direction | Sample-overlap mis-specification OR sign error in munging | Re-munge with same allele convention; re-run `ldsc()` |
| MTAG MaxFDR > 5%, GenomicSEM with Q_SNP works | MTAG assumption violated | Prefer GenomicSEM as primary |
| One-trait GWAS sig but common-factor not | Trait-specific architecture | Don't force into common-factor frame |

**Operational rule for publication:** A common-factor SNP claim requires (1) factor p < 5e-8, (2) Q_SNP p > 0.05 / N_factor_SNPs (non-heterogeneous), and (3) replication in an independent set of traits or cohorts. Trait-specific SNPs from MTAG require MaxFDR < 5% for the trait. Reporting only the factor effect without Q_SNP is the most common reviewer-flagged error.

## MTAG vs GenomicSEM Common-Factor GWAS

Both methods exploit genetic correlation among input GWAS, but their goals and outputs differ.

| Property | MTAG | GenomicSEM commonfactorGWAS |
|----------|------|------------------------------|
| Output | Per-trait shrunk z-scores | SNP effect on latent factor |
| Sample-overlap handling | Bivariate LDSC intercept | Full LDSC sampling-covariance matrix V |
| Heterogeneity diagnostic | MaxFDR (Turley 2018) | Q_SNP (Grotzinger 2019) |
| Interpretation | "Boosted power for trait k" | "Effect on what the traits share" |
| Min traits | 2 | 3 (otherwise factor not identified) |
| Best when | Power-boost an individual trait | Common factor hypothesized |

Both depend on accurate sampling covariance. MTAG fails (MaxFDR > 5%) under the same heterogeneity that produces large Q_SNP in GenomicSEM. The two methods should be reported together when the prior on a common factor is non-trivial; agreement increases confidence, disagreement points to architecture-specific SNPs. Prefer `commonfactorGWAS` over MTAG when the traits fit a common factor (CFI >= 0.95): it models heterogeneity explicitly through Q_SNP.

## Per-Method Failure Modes

### MTAG MaxFDR > 5%

**Trigger:** Running MTAG on traits with low pairwise genetic correlation or with one trait that has a very different architecture.

**Mechanism:** MTAG assumes a homogeneous variance-covariance structure across SNPs. When heterogeneity dominates, the empirical-Bayes shrinkage can over-claim SNPs in the focal trait. Turley 2018 defines MaxFDR as the maximum estimated false discovery rate under worst-case heterogeneity; > 5% invalidates the published trait-specific summary statistics.

**Symptom:** MTAG output file reports `maxFDR` > 0.05; per-trait MTAG hits don't replicate in independent cohorts.

**Fix:** Check pairwise rg via LDSC; if any pair is < 0.7, MTAG is risky. Drop the most heterogeneous trait and re-run. Alternatively, switch to GenomicSEM's `commonfactorGWAS` which models heterogeneity explicitly via Q_SNP.
