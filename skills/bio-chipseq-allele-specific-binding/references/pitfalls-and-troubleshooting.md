# ASB pitfalls and troubleshooting

Use this reference to diagnose unexpected allelic ratios or tool failures. Preserve excluded loci and diagnostics so a failed assumption remains visible.

## Universal confounders

### Reference-allele mapping bias

**Trigger:** BaalChIP, a beta-binomial test, or a chi-squared test is run directly on ordinary reference-aligned reads.

**Mechanism:** Reads carrying the reference allele often align more readily than reads carrying an alternate allele, creating false REF-favoring calls.

**Symptoms:** Significant sites cluster toward REF; the tested-site ratio distribution is systematically above 0.5.

**Action:** Apply WASP with sample-specific variant inputs, use RASQUAL's modeled bias path, or use a validated personalized-genome workflow. Recount and re-evaluate attrition.

### Imprinted loci

**Trigger:** Strong ASB is concentrated at H19, IGF2, MEG3, MEG8, KCNQ1OT1, or another known imprinted region.

**Mechanism:** Parent-of-origin expression or chromatin state produces genuine constitutive allele skew that is not evidence of a local binding variant.

**Action:** Filter or separately label imprinted loci before ordinary ASB interpretation. Record the catalog and genome-build mapping used.

### X-inactivation

**Trigger:** Female samples show widespread chrX imbalance.

**Mechanism:** X-inactivation creates allele skew across X-linked loci.

**Action:** Exclude chrX from the autosomal analysis or use a separate, explicitly X-aware interpretation.

### Copy-number imbalance

**Trigger:** Cancer-sample allelic ratios are bimodal or significant calls cluster in altered segments.

**Mechanism:** Allele-specific gains and losses change effective allele dose, so raw count imbalance mixes copy number with binding.

**Action:** Use BaalChIP with measured per-variant RAF or matching gDNA, or exclude altered segments. Do not treat ordinary population AF or an unconsumed CN-segment BED as correction. Confirm that every input's assembly and contig style match.

## Tool-specific failure modes

### WASP: variant-panel mismatch

**Trigger:** The WASP SNP file comes from a population panel rather than the sample.

**Mechanism:** Reads overlapping sample-specific variants absent from the panel are never swapped and remapped.

**Symptom:** Residual reference bias at sample-specific heterozygous sites.

**Action:** Rebuild WASP inputs from the sample's genotype VCF, or choose a validated alternative that handles those sites.

### WASP: high read loss

**Trigger:** More than the expected project-specific fraction of reads is discarded, especially above 40%.

**Mechanism:** Reads spanning several heterozygous SNPs must remap consistently after allele swaps; repetitive or variant-dense regions can fail disproportionately.

**Action:** Verify paired-end handling, SNP-table construction, reference build, intermediate counts, and remap settings. If the attrition is real, report the reduced power; consider a personalized-genome or modeled-bias strategy rather than bypassing correction.

### RASQUAL: convergence failure

**Trigger:** Features have few informative heterozygous SNPs, poor imputation, or strong copy-number imbalance.

**Mechanism:** Sparse allelic information underspecifies the feature model.

**Action:** Require well-imputed feature SNPs, combine justified replicates, narrow and verify feature windows, or use an analysis suited to sparse cancer data.

### BaalChIP: no overlap with peaks

**Cause:** Heterozygous variants and peaks use different assemblies, coordinate conventions, or chromosome prefixes.

**Action:** Compare BAM contigs and representative intervals from every input; standardize only after recording the transformation.

### BaalChIP: copy-number resource mismatch

**Trigger:** Copy-number intervals use `X` while BAM/variant inputs use `chrX`, or otherwise fail to overlap.

**Symptom:** Calls in altered regions remain strongly bimodal.

**Action:** Verify exact overlap for exclusion/stratification, but separately confirm that correction enters BaalChIP through a measured `RAF` column or group-named `CorrectWithgDNA` BAM. The tested 1.36.0 API does not consume a detached CNV BED.

### AlleleSeq: insufficient phasing

**Trigger:** An unphased or poorly phased VCF is used to create the diploid genome.

**Mechanism:** Haplotype alleles are assigned arbitrarily, destroying maternal/paternal interpretation and fragment consistency.

**Action:** Phase with trio or read-backed methods such as HapCUT2 or WhatsHap, and report switch-error limitations.

### AlleleSeq: diploid genome too large

**Cause:** The input VCF includes many or large structural variants.

**Action:** Use an explicitly small-variant-only VCF or exclude SV-rich regions, while recording the resulting scope.

## Error and symptom table

| Error or symptom | Likely cause | Next check |
|---|---|---|
| WASP `find_intersecting_snps.py` fails | SNP HDF5/text layout does not match the release | Rebuild with that release's `snp2h5` or extraction utility |
| Calls remain REF-skewed | WASP omitted, sample-specific variants absent, or remap path incomplete | Inspect retained-read ratios and input SNP coverage |
| All calls occur at imprinted loci | Imprinting filter omitted | Apply and document an assembly-matched imprinted-locus set |
| Many calls occur on female chrX | X-inactivation | Separate chrX from the autosomal analysis |
| Cancer ratios are bimodal | Allele-specific copy number ignored | Overlay calls with CN segments and verify model ingestion |
| RASQUAL runs out of memory | Cis-window or feature batch too broad | Narrow the tabix region and `-l`/`-m`; chunk features |
| BaalChIP reports no hetSNPs in peaks | Assembly/prefix/coordinate mismatch | Compare exact intervals across BAM, variant, and peak inputs |
| Significant table is empty | Low depth, aggressive attrition, incompatible inputs, or conservative model | Trace counts per stage before changing thresholds |
| ASB disagrees with chromBPNet | Low sample coverage, cell-state mismatch, or model extrapolation | Inspect read depth and model ensemble/context |

## Reconciliation patterns

| Pattern | Interpretation to test | Action |
|---|---|---|
| WASP applied but REF skew remains | Sample variants were absent from WASP inputs or remap was incomplete | Audit SNP-table coverage and intermediate files |
| BaalChIP and RASQUAL disagree at sparse sites | Models respond differently to sparse evidence | Inspect counts, posterior/uncertainty, `phi`, and CN context |
| ASB at an imprinted locus | Real biology, but not ordinary differential binding evidence | Label separately or exclude from the primary result |
| ASB on female chrX | X-inactivation may dominate | Use a separate X-aware analysis |
| ASB inside a CN-altered region | Allele dose may dominate | Verify measured RAF/gDNA correction or exclude the segment |
| Strong predicted effect without ASB | The prediction is counterfactual; the sample may lack depth or the relevant state | Report both with coverage and context |

## Power and reporting cautions

- Mapping-bias correction can materially reduce read depth; measure actual loss rather than promising a fixed percentage.
- Determine depth and power from the selected model, dispersion, test family, replicate structure, and expected effect; do not promote an unverified fixed read threshold to a universal requirement.
- Combining replicates can increase counts but can also hide replicate heterogeneity. Preserve replicate-level diagnostics.
- Do not relax filters merely to recover a desired number of hits. Trace whether losses arise from alignment bias, peak restriction, depth, blacklist, imprinting, chrX, CN, or model uncertainty.
