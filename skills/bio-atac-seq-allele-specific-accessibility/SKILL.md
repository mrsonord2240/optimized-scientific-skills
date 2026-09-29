---
name: bio-atac-seq-allele-specific-accessibility
description: Detect allele-specific chromatin accessibility from ATAC-seq using WASP, GATK ASEReadCounter, or RASQUAL. Use when mapping cis-regulatory genetic variants from heterozygous SNPs, validating GWAS variant function with allelic imbalance, building a phased within-peak ASE analysis, or detecting reference-allele mapping bias before downstream analysis.
tool_type: mixed
primary_tool: WASP
license: MIT
author: GPTomics
---

# Allele-Specific Accessibility

Determine whether a heterozygous SNP changes chromatin accessibility on its
allele by counting reference- and alternative-supporting ATAC reads within an
individual. Treat this as a binomial allelic-imbalance problem, not as ordinary
cohort differential accessibility. In the source workflow, per-site power is
driven by allelic read depth and informative heterozygous sites rather than by
treating the number of individuals as the binomial sample size.

## Verify the environment

Reference examples target WASP 0.3.4+, GATK 4.4+, RASQUAL 1.1+, samtools
1.19+, bcftools 1.19+, vcftools 0.1.16+, plink 2.00+, bowtie2 2.5+,
scipy 1.11+, pandas 2+, and pybedtools 0.10+. MatrixEQTL and QuASAR are
method-selection routes only in this Skill; no executable candidate workflow
for either is shipped or claimed.

Before using a command or code pattern, inspect the installed interface:

- CLI: run `<tool> --version` and `<tool> --help`.
- Python: run `pip show <package>` and inspect the called function signature.
- R: run `packageVersion('<pkg>')` and open the called function's help page.

If the installed interface differs, adapt the invocation rather than retrying
unchanged. Read [setup and inputs](references/setup-and-inputs.md) for
dependencies, required files, and preflight checks.

## Choose the analysis

Use the smallest workflow that answers the question:

| Setting | Pipeline |
|---|---|
| One individual with genotypes | WASP + GATK ASEReadCounter, then per-SNP and per-peak ASE |
| Cohort of 20-100 | WASP + RASQUAL for joint total-count and allelic evidence |
| Cohort of at least 100 | Route to a separately specified and validated cohort peak-count model; this Skill does not ship one |
| No genotypes | Route to a separately specified and validated genotype/ASE inference workflow; this Skill does not ship one |
| GWAS-variant validation | Identify heterozygous carriers and aggregate ASE at the variant |
| Trio or quartet | Phase within the family, then analyze each individual |

Read [method selection and failure modes](references/method-selection-and-failures.md)
before choosing among these branches. For RASQUAL, also read
[RASQUAL joint modeling](references/rasqual-joint-modeling.md).

## Run the core workflow

1. Confirm that the ATAC BAM is deduplicated and mapping-quality filtered, the
   VCF belongs to the same sample, the reference FASTA matches the aligner
   index, and the peak set uses the same assembly and chromosome naming.
2. Resolve exactly one BAM `@RG/SM` sample, require the same sample in the VCF,
   subset to it, and only then restrict to biallelic heterozygous SNPs. Phasing
   and a non-missing phase-set identifier are required for within-peak pooling;
   per-site GATK counting itself does not require phase.
3. Build WASP SNP and haplotype HDF5 tables with `snp2h5`.
4. Use WASP to identify reads overlapping heterozygous SNPs, swap alleles,
   realign the swapped reads, and retain only reads that map identically.
5. Merge the retained remapped reads with reads that did not require remapping,
   then sort the merged BAM, add a validated sample read group to orphan reads,
   and index it. Do not count alleles from the original, uncorrected BAM.
6. Run GATK ASEReadCounter on the WASP-filtered BAM and filtered VCF. Stop if
   GATK writes no data rows even when its process exits successfully.
7. Analyze per-SNP counts or aggregate SNPs only after orienting every count to
   the same phased haplotype. Keep different phase sets separate. Apply a
   binomial test against 0.5, correct p-values with BH FDR, and report both
   haplotype counts, coverage, effect size, SNP count, and call status.
8. For a moderate cohort, use the validated per-feature RASQUAL driver and
   apply its BH summary across every reported feature-variant test.
   Large-cohort and no-genotype methods require a separate executable analysis
   contract.

The complete runnable single-sample workflow is
[`scripts/wasp_ase_pipeline.sh`](scripts/wasp_ase_pipeline.sh). It invokes
[`scripts/aggregate_peak_ase.py`](scripts/aggregate_peak_ase.py) for the final
within-peak aggregation step. The bounded RASQUAL driver is
[`scripts/run_rasqual_features.py`](scripts/run_rasqual_features.py).

## Enforce the mapping-bias invariant

Reference alleles align more easily to a reference genome than alternative
alleles. This can shift null reference fractions above 0.5 and create false
allelic imbalance. Always run WASP before GATK ASEReadCounter or any equivalent
allele counter; there are no exceptions in this workflow.

Verify the correction empirically by plotting the reference-allele fraction
across presumed-null sites. A systematic reference excess means the counts are
not safe for inference.

## Interpret within-peak results

Use these source defaults unless the study protocol specifies otherwise:

- Require `totalCount >= 30` before per-SNP interpretation.
- Require at least two heterozygous SNPs for a pooled peak call.
- Call statistical significance at BH-adjusted `p < 0.05`.
- Require `abs(ref_frac - 0.5) >= 0.2` for the source's biologically meaningful
  30:70 imbalance threshold.
- Treat lower-depth or smaller-effect results as underpowered rather than as
  evidence of balance.

## Report evidence and limitations

For every result, record software versions, reference assembly, genotype
phasing method, WASP filtering, read/base-quality thresholds, coverage filters,
multiple-testing correction, and effect-size threshold. Distinguish per-SNP,
per-peak, per-individual, and cohort-level evidence.

Report WASP-filtered ASE with its adjusted p-value and effect size. For cohort
evidence, define the complete feature/variant test family in advance and use a
study-specific FDR or permutation procedure; a fixed raw p-value is not a
substitute for multiplicity control. MPRA or CRISPRi-FlowFISH evidence is a
separate validation layer and is not produced by this Skill.

Consult [method selection and failure modes](references/method-selection-and-failures.md)
when tools disagree, coverage is sparse, RASQUAL returns missing results, or a
claim fails to replicate across individuals.

## Resource map

- [Setup and inputs](references/setup-and-inputs.md): installation clues,
  required inputs, compatibility checks, and sample-identity preflight.
- [Method selection and failure modes](references/method-selection-and-failures.md):
  tool taxonomy, thresholds, failure diagnosis, reconciliation, literature,
  and related Skills.
- [RASQUAL joint modeling](references/rasqual-joint-modeling.md): binary matrix
  preparation, streamed-VCF invocation, and flag meanings.
- [Provenance](references/provenance.md): pinned provider source, authorship,
  license, and source-file checksums.

