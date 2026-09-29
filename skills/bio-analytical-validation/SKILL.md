---
name: bio-analytical-validation
description: Treats a ctDNA assay as a molecule-counting experiment at the Poisson edge and builds its analytical-validation case the measurement-science way. Covers the genome-equivalent currency (~303 haploid copies/ng under a 3.3 pg convention), the lambda = input_GE x VAF sampling ceiling (lambda~=2.996 for 95% template presence), the error-suppression ladder (raw NGS ~1e-3 -> single-strand UMI ~1e-4/1e-5 -> duplex <1e-7), the CLSI EP17 LoB/LoD/LoD95/LoQ framework, the distinction between theoretical multi-locus sampling bounds and empirically validated panel LoD95, contrived/SEQC2 reference standards, and honest LoD reporting conditioned on input mass + consensus depth + replicate detection rate. Use when stating or trusting a sensitivity claim, designing a dilution-series validation, deciding how many genome equivalents are needed at a target VAF, distinguishing sampling calculations from assay LoD, or auditing a "detects 0.1% VAF" claim.
tool_type: python
primary_tool: scipy
---

## Version Compatibility

Reference examples tested with: numpy 1.26+, scipy 1.12+, statsmodels 0.14+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Analytical Validation and Detection Limits

**"What is the real limit of detection of my ctDNA assay, and can I trust the number I am about to report?"** -> Quantify the Poisson sampling ceiling, the error-suppression floor, and the LoB/LoD/LoQ that together define a defensible sensitivity claim.
- Python: `scipy.stats.poisson` for detection-probability math, `scipy.stats.norm` for CLSI LoB/LoD, `statsmodels` Probit/Logit for a dilution-series LoD95 fit.

## The Single Most Important Modern Insight -- LoD Is Set by Genome Equivalents Sampled and Error Suppression, NOT by Sequencing Depth or the Caller

A ctDNA assay is a molecule-counting experiment at the Poisson edge. The mutant signal is a fixed, tiny number of physical template molecules in the tube, and the limit of detection is governed by two ceilings: how many genome equivalents were sampled (Poisson), and how low the background error floor was driven (error suppression). Under the explicit 3.3 pg per haploid-genome convention, 1 ng contains about 303 haploid genome equivalents (1,000 pg / 3.3 pg); laboratories using another conversion must state it. The expected mutant-molecule count is lambda = input_GE x VAF. A 0.1% variant on 1,000 GE (~3.3 ng) has lambda = 1, so e^-1 ~= 37% of the time the mutant template was never in the tube and a perfect sequencer detects nothing. Past the point where every input molecule has been read once (sampling saturation, visible as a deduplication plateau in UMI families), additional read depth re-sequences the same physical molecules and adds zero information. Reporting an LoD as a bare VAF -- with no input mass, no unique-molecule (consensus) depth, no replicate detection rate -- is reporting an undefined quantity.

The second ceiling is the per-base background error rate, which sets the VAF floor independently: a 0.1% variant cannot be distinguished from noise if the assay manufactures that base at 0.1%. Error suppression is a ladder (raw NGS ~1e-3 -> single-strand UMI consensus ~1e-4/1e-5 -> duplex <1e-7), and single-strand consensus does NOT remove template-resident damage (C->T deamination, G->T 8-oxoG) because every PCR copy of that strand inherits the lesion -- only duplex strand-concordance catches it. The achieved LoD is the *worse* of the sampling and error ceilings. Multi-locus integration can improve assay sensitivity, but a binomial calculation over ideal Poisson template presence is only an optimistic theoretical sampling lower bound. An achieved panel LoD95 requires full-workflow replicate data that include recovery, consensus depth, locus-specific background, false positives, and the actual calling rule.

## Methods Landscape

| Concept | Definition | Source |
|---------|------------|--------|
| LoB (Limit of Blank) | Highest signal expected from an analyte-free blank (95th pct): LoB = mean_blank + 1.645*SD_blank; the false-positive anchor on true negatives | CLSI EP17-A2 |
| LoD (Limit of Detection) | Lowest level reliably distinguishable from LoB: LoD = LoB + 1.645*SD_low; a sample at LoD is detected ~95% of the time | CLSI EP17-A2 |
| LoD95 | The concentration/VAF where detection probability = 95%; a point on a probit/logistic detection curve, not a separate definition | CLSI EP17-A2; Newman 2016 |
| LoQ (Limit of Quantitation) | Lowest level measurable with stated precision (e.g. CV<=20%); LoQ >= LoD always, so a "VAF" near the floor is detectable but not trustworthy | CLSI EP17-A2 |
| Per-locus LoD | Single-variant LoD; sampling- and error-limited (~0.05-0.1% VAF typical) | Newman 2014/2016 |
| Theoretical multi-locus sampling bound | Binomial over independent, equal-VAF Poisson template-presence probabilities under a >=k-of-N rule; an optimistic lower bound, not assay LoD95 | Analytical model; assumptions must accompany every result |
| Empirical panel LoD95 | VAF where the complete assay and fixed panel calling rule detects >=95% of full-workflow replicates, reported with uncertainty | CLSI EP17-style detection study; panel rule fixed before fitting |
| HCC1395 truth set | Tumor HCC1395 (ATCC CRL-2324) and matched normal HCC1395BL (ATCC CRL-2325); SEQC2 somatic call set v1.2 and SRA SRP162370 | Fang 2021; DOI 10.1038/s41587-021-00993-6 |
| SEQC2 oncopanel/ctDNA materials | Sample A is an equal-mass pool of ten UHRR cancer-cell-line DNAs; Sample B is Agilent male reference DNA 5190-8848; A/B dilutions Df/Ef are fragmented to 130-170 bp | Jones 2021, DOI 10.1186/s13059-021-02316-z; Deveson 2021, DOI 10.1038/s41587-021-00857-z |

## Decision Tree by Scenario

| Scenario | Recommended | Why |
|----------|-------------|-----|
| "How many GE for 95% template presence at VAF X?" | Solve lambda = input_GE x VAF = -ln(0.05), so input_GE ~= 2.996/VAF | This is sampling only; ~29,957 GE (~98.9 ng at 303.03 GE/ng) for one 1e-4 variant, before recovery or calling losses |
| "Why is more depth not helping?" | Report unique (consensus) molecular coverage, not raw depth; check the dedup plateau | Past sampling saturation, depth re-reads the same molecules; the ceiling is GE in the tube |
| Single hotspot vs bespoke panel for low VAF | Use the multi-locus calculation only to explore an ideal sampling bound; validate the chosen panel and >=k-of-N rule empirically | Real loci differ in recovery and background, and detections are not guaranteed independent |
| Reporting an LoD | Condition on input mass (GE) + consensus depth + replicate detection rate (e.g. "LoD95 0.1% VAF at 30 ng / 2x duplex / 95% of 20 replicates") | A bare VAF omits the input mass, the unique depth, and per-locus vs integrated -- it is undefined |
| Estimating LoD95 from a dilution series | Use replicated levels spanning below and above 95%; reject separation/nonpositive slope, and report inverse-prediction CI plus convergence and slope diagnostics | Binary probit/logistic models require identified, bracketed data; a point estimate alone is insufficient |
| Distinguishing detection from quantitation | Set LoD for yes/no calls; set LoQ (CV<=20%) separately for any reported VAF/TF | MRD calls are binary and can sit far below LoQ; a near-floor VAF number is not quantitative |
| Validating somatic callers | Use HCC1395 (ATCC CRL-2324) with HCC1395BL (CRL-2325) and the SEQC2 v1.2 call set | This paired tumor-normal truth-set program is distinct from Sample A materials |
| Validating oncopanel/ctDNA workflows | Use the appropriate SEQC2 Sample A/B-derived dilution (for example Df/Ef at 130-170 bp) and record its exact study material | Sample A is a ten-cell-line pool, not HCC1395; commutability with patient plasma remains a caveat |

## Genome-Equivalent and Poisson Detection Calculator

**Goal:** Convert an input mass and target VAF into an expected mutant-molecule count and a detection probability, so a sensitivity claim is anchored to molecules rather than to a VAF alone.

**Approach:** Under the default 3.3 pg convention, convert ng to 303.03 haploid genome equivalents/ng, set lambda = input_GE x VAF, and read template-presence probability as a Poisson tail P(X >= k). Invert the continuous Poisson mean for the requested probability. Label every result sampling-only: it excludes recovery, consensus depth, background error, false positives, and calling behavior. Override `ge_per_ng` only with a documented laboratory convention.

See [scripts/](scripts/) for implementation:
- `ge_and_poisson.py` - genome equivalents, Poisson detection probability, and sampling-detection thresholds

## LoB and LoD95 from a Dilution Series (CLSI EP17 style)

**Goal:** Estimate the VAF at which the assay detects 95% of the time, from a contrived dilution series, and anchor it to the blank-derived false-positive floor.

**Approach:** Compute the Gaussian-shortcut LoB from at least two finite blank replicates (mean + 1.645*SD), then fit a probit GLM of binary detection on log10(VAF). Require at least three replicated VAF levels, mixed binary outcomes, and empirical detection rates that bracket 0.95. Reject separated, non-converged, nonpositive-slope, and physically impossible fits. Report the LoD95 inverse-prediction confidence interval and diagnostics returned by `lod95_probit`; a point estimate without them is not acceptable.

See [scripts/](scripts/) for implementation:
- `lod95_probit.py` - limit of blank and LoD95 probit fit

## Theoretical Multi-Locus Sampling Lower Bound

**Goal:** Explore the best-case sampling boundary implied by input mass and a >=k-of-N rule without claiming achieved assay sensitivity.

**Approach:** Treat each tracked locus as an independent Poisson sampler at the same VAF and assume perfect recovery and detection, no consensus-depth failure, no background errors or false positives, and no locus-to-locus heterogeneity. A panel-positive call requires at least k templates present. The result is an optimistic VAF lower bound. Do not call it LoD95; validate the complete workflow on replicated materials to establish panel LoD95.

See [scripts/](scripts/) for implementation:
- `panel_integrated_lod.py` - structured sampling-only probability and theoretical VAF lower-bound calculations with assumptions

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| 303.03 haploid genome equivalents per ng cfDNA | Explicit 3.3 pg per haploid-genome convention | 1,000 pg / 3.3 pg = 303.03. Other conventions are permitted only when stated and passed as `ge_per_ng`; 330 does not follow from 3.3 pg or 6.6 pg without an additional approximation. |
| lambda = input_GE x VAF; lambda >= 3 for ~95% sampling-detection | Poisson, 1 - e^-3 = 0.95 | Below lambda~3 the mutant template is often simply absent from the tube regardless of sequencing |
| Raw NGS error floor ~1e-3 | Schmitt 2012 context; field consensus | Sets the per-base VAF floor before any consensus; a global VAF cutoff above this is noise-limited |
| Single-strand UMI consensus ~1e-4 to 1e-5 | Newman 2014/2016 (CAPP-Seq/iDES) | Majority-vote within a UMI family erases PCR/sequencing error not shared across the family |
| Duplex sequencing <1e-7 (theory <1/1e9 nt) | Schmitt 2012 *PNAS* 109:14508 | Requires the variant on BOTH original strands; independent strand errors cannot agree |
| iDES adds ~3-15x over baseline; ctDNA to ~4e-5 | Newman 2016 *Nat Biotechnol* 34:547 | Position/trinucleotide background model subtracts stereotyped artifacts per locus |
| ichorCNA tumor-fraction floor ~3% | Adalsteinsson 2017 *Nat Commun* 8:1324 | Copy-number-based TF estimation; sWGS/ULP-WGS cannot resolve TF below ~3% -- an LoD, not a VAF |
| Multi-locus integration can improve sensitivity | Reinert 2019 *JAMA Oncol* 5:1124 | Clinical performance is empirical; the shipped independent-Poisson calculation is only a sampling lower bound and cannot establish ppm LoD95 |
| LoQ >= LoD (e.g. CV<=20% for quantitation) | CLSI EP17-A2 | Detection (binary) is easier than quantitation (continuous); near-floor VAFs are not trustworthy numbers |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| "Assay detects 0.1% VAF" with no input mass | VAF reported as a standalone sensitivity spec | Condition the LoD on input GE + consensus depth + replicate detection rate; 0.1% on 100 GE is noise |
| Buying more sequencing depth to improve sensitivity | Conflating read depth with molecule count | Past the dedup plateau the assay is sampling-saturated; add plasma volume / conversion efficiency, not depth |
| Sampling-only k-of-N result quoted as panel LoD95 | Treating template presence as complete assay detection | Call it a theoretical sampling VAF lower bound and state all assumptions; use replicated full-workflow detections for assay LoD95 |
| VAF used as the sensitivity unit | Omitting the molecule count behind the fraction | Pair every VAF with input GE; lambda = GE x VAF is the quantity that determines detection |
| Single-strand UMI assumed to remove damage artifacts | Template-resident C->T/G->T inherited by every copy | Use duplex strand-concordance for sub-1e-5 claims; single-strand votes unanimously for the lesion |
| Reporting a near-floor VAF as a measured value | Confusing LoD (detect) with LoQ (quantify) | Quantitative VAF/TF only at/above LoQ (CV<=20%); below it report detected/not-detected |
| Global VAF cutoff across all loci | Background error is position/context-dependent | Use a per-locus background model (iDES-style); a flat threshold loses sensitivity and specificity |

## References

- Diehl F, Schmidt K, Choti MA, et al. 2008. Circulating mutant DNA to assess tumor dynamics. *Nat Med* 14:985-990. -- ctDNA half-life ~114 min; molecule-counting framing of tumor dynamics.
- Schmitt MW, Kennedy SR, Salk JJ, et al. 2012. Detection of ultra-rare mutations by next-generation sequencing. *PNAS* 109:14508-14513. -- Duplex sequencing; theoretical error floor <1 per 1e9 nt.
- Newman AM, Bratman SV, To J, et al. 2014. An ultrasensitive method for quantitating circulating tumor DNA with broad patient coverage. *Nat Med* 20:548-554. -- CAPP-Seq; UMI-consensus error suppression.
- Newman AM, Lovejoy AF, Klass DM, et al. 2016. Integrated digital error suppression for improved detection of circulating tumor DNA. *Nat Biotechnol* 34:547-555. -- iDES; ~3-15x gain; ctDNA to ~4e-5.
- Razavi P, Li BT, Brown DN, et al. 2019. High-intensity sequencing reveals the sources of plasma circulating cell-free DNA variants. *Nat Med* 25:1928-1937. -- CHIP as the dominant non-tumor signal in the LoB blank.
- Adalsteinsson VA, Ha G, Freeman SS, et al. 2017. Scalable whole-exome sequencing of cell-free DNA reveals high concordance with metastatic tumors. *Nat Commun* 8:1324. -- ichorCNA; copy-number tumor-fraction floor ~3%.
- Reinert T, Henriksen TV, Christensen E, et al. 2019. Analysis of plasma cell-free DNA by ultradeep sequencing in patients with stages I to III colorectal cancer. *JAMA Oncol* 5:1124-1131. -- Signatera; 16-variant integration; >=2-of-N positivity.
- Fang LT, Zhu B, Zhao Y, et al.; SEQC2 Consortium. 2021. Establishing community reference samples, data and call sets for benchmarking cancer mutation detection using whole-genome sequencing. *Nat Biotechnol* 39:1151-1160. DOI 10.1038/s41587-021-00993-6. -- HCC1395 (ATCC CRL-2324), HCC1395BL (CRL-2325), call set v1.2, SRA SRP162370; not the SEQC2 Sample A pool.
- Jones W, et al.; SEQC2 Consortium. 2021. A verified genomic reference sample for assessing performance of cancer panels detecting small variants of low allele frequency. *Genome Biol* 22:111. DOI 10.1186/s13059-021-02316-z. -- Sample A ten-cell-line pool and Sample B (Agilent 5190-8848).
- SEQC2 Oncopanel Sequencing Working Group. 2021. Evaluating the analytical validity of circulating tumor DNA sequencing assays for precision oncology. *Nat Biotechnol* 39:1115-1128. DOI 10.1038/s41587-021-00857-z. -- Sample A/B-derived liquid-biopsy materials and fragmented Df/Ef/Ff; BioProject PRJNA677999.
- CLSI EP17-A2. 2012. Evaluation of Detection Capability for Clinical Laboratory Measurement Procedures; Approved Guideline -- Second Edition. Clinical and Laboratory Standards Institute. -- Governing LoB/LoD/LoQ definitions.

## Related Skills

- ctdna-mutation-detection - applies these limits to low-VAF somatic calls
- longitudinal-monitoring - per-timepoint LoD and left-censoring of undetectable samples
- tumor-fraction-estimation - the ~3% CNA-based detection floor as an LoD
- experimental-design/multiple-testing - repeated-surveillance specificity and FDR
- clinical-biostatistics/power-and-sample-size - validation-study design
