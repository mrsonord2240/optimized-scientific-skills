# Background

Not needed to run an analysis. This explains choices the routes and scripts already make.

## Governing principle

Filtering decides which errors a callset keeps, not whether they exist. The two site-level paradigms fail in opposite regimes: VQSR (a ratio of two learned Gaussian-mixture densities over annotation space) collapses on small or exome cohorts; static hard thresholds discard real variants at scale. Site-level filtering ("is this site real?") and genotype-level filtering ("is this sample's genotype trustworthy?") are orthogonal, and both are needed, site first. SNPs and indels have different error processes (base-calling and strand against alignment ambiguity in repeats) and different truth resources, so they are always filtered separately and then merged. None of these mistakes throws an error: the VCF stays valid while the numbers are silently wrong.

## VQSR mechanism

Fit a Gaussian mixture to the annotation profile of known-true sites (positive model), bootstrap a second mixture on the low-probability-tail artifact sites (negative model), and score each variant by VQSLOD = log(P(annotations | positive) / P(annotations | negative)). A "99.7 tranche" is the VQSLOD cutoff that retains 99.7% of the truth-set sites. VQSR estimates full covariance matrices in 6-8 dimensions, so it needs tens of thousands of variants: as practical GATK convention a single deep WGS or about 30 or more jointly genotyped exomes. Indels are about 10 times rarer than SNPs, so their model fails first.

## Quality metric rationale

| Metric | Why |
|--------|-----|
| QD | QUAL divided by alt-supporting depth. Raw QUAL grows with coverage, so a 500x artifact can post a large QUAL; QD removes that. Bimodal in practice: real variants about 12-35, artifacts near 0. The workhorse, not QUAL |
| FS | Phred-scaled Fisher's exact p-value for strand bias. Real variants are strand-symmetric, many artifacts are not. Breaks down at exon and read ends, where SOR takes over |
| SOR | Symmetric odds ratio; tolerates the legitimate strand imbalance at exon and read ends where FS gives false positives. Complements FS |
| MQ | RMS mapping quality at the site. Low means ambiguous mapping (repeats, paralogs, segmental duplications) |
| MQRankSum | Rank-sum of mapping quality, alt against ref reads. Strongly negative means alt reads map worse. Missing at hom-alt sites |
| ReadPosRankSum | Rank-sum of within-read position. Strongly negative means alt clusters at read ends. Missing at hom-alt sites |
| DP | Extreme depth (above 2x or below 0.3x the mean) suggests collapsed repeats or poor capture; never a VQSR annotation on exomes |
| GQ | Genotype level. Phred confidence in the called genotype; GQ 20 is 99% |

The "missing passes" rule: the rank-sum tests compare ref-supporting and alt-supporting reads, so they are defined only at heterozygous sites. GATK `VariantFiltration` fires a filter only when the value is present and violates the cutoff.

## Other pitfalls

- Re-applying GATK hard filters to DeepVariant or DRAGEN output: the annotations may be absent or differently distributed.
- Choosing thresholds without looking at the annotation distributions.
- Ti/Tv rationale: random errors have Ti/Tv about 0.5, diluting the signal; WES is higher than WGS because of CpG-transition-rich coding enrichment.
- Genotype filtering before cohort metrics: HWE on low-GQ genotypes manufactures deviation.

## Open points for audit

- Source cyvcf2 snippet included a DP floor (10) its own text said "mirrors the bcftools expression term for term", but the bcftools expression has no DP term. Kept as `min_dp=10` option in `scripts/cyvcf2_filter.py`.
- Source cyvcf2 snippet applied SNP thresholds to every variant (its own listed pitfall); the script now splits by type.
- Source somatic expression did not require `FILTER=PASS`; kept as the default, `pass_only=1` added.
- GATK 4.6.1.0 `gatkdoc` lists `ScoreVariantAnnotations` as beta; it has no page for `CNNScoreVariants` or `NVScoreVariants`, so the deprecation claim is carried from the source and unconfirmed here.
