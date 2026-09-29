# Method selection and failure modes

Use this reference when choosing a branch, interpreting sparse evidence, or
reconciling outputs from multiple tools.

## Tool taxonomy

| Tool | Method and input | Strength | Main limitation |
|---|---|---|---|
| WASP | Allele-swap remapping from BAM, VCF, and reference | Controls reference-allele mapping bias and is mandatory before allele counting | Realignment is slow on deep data |
| GATK ASEReadCounter | REF/ALT counts at known heterozygous sites | Mature single-sample counter in the GATK ecosystem | Does not correct mapping bias; run only after WASP |
| RASQUAL | Joint total-count and allelic model from peak counts and a streamed VCF | Higher power for moderate cohorts and models genotype uncertainty | Non-standard setup and per-feature execution |
| QuASAR | External route: joint genotype and ASE inference from sequencing reads | Can support studies without external genotypes | This Skill ships no executable QuASAR data/model/QC/output contract |
| MatrixEQTL | External route: linear association between genotypes and peak counts | Standard cohort-level caQTL analysis | This Skill ships no executable MatrixEQTL matrices/covariates/model/output contract |
| Bayesian beta-binomial models | Overdispersed allelic-count model | Represents overdispersion explicitly | Less standardized in this workflow |

Modern caQTL studies commonly combine WASP with RASQUAL at moderate sample
sizes, or WASP plus GATK and a linear caQTL model in large cohorts. Verify the
current primary literature and installed implementations before freezing a
pipeline.

## Decision details

- For one genotyped individual, count alleles with GATK after WASP and aggregate
  compatible alleles within peaks to improve power.
- For 20-100 individuals, consider RASQUAL so total accessibility and allelic
  imbalance contribute to one joint test.
- For at least 100 individuals, route to a separately specified cohort
  peak-count analysis with explicit matrices, covariates, QC, cis windows, and
  multiplicity control; retain per-individual ASE as independent cis evidence.
- Without genotypes, route to a separately specified genotype/ASE inference
  workflow with an executable read-to-count contract and fit diagnostics. This
  Skill does not provide a runnable QuASAR branch.
- For a GWAS candidate, identify heterozygous carriers, assess the variant in
  each carrier, and meta-analyze rather than treating reads from different
  people as one binomial sample.
- For families, phase before haplotype-oriented pooling and retain the
  individual as the unit of analysis.
- For an iPSC line with known SNPs, use the single-individual branch to test
  allele imbalance and, when phase is known, parental origin.

## Power and thresholds

Per-SNP ATAC coverage is often only 10-100 reads. Ten reads provide almost no
power for a moderate imbalance, and at 30 reads only large shifts are readily
detectable. Pool replicates when scientifically valid, aggregate compatible
SNPs within a peak, or combine ASE with cohort caQTL evidence.

The provider defaults are:

| Quantity | Default |
|---|---|
| Per-SNP coverage | `totalCount >= 30` |
| SNPs per pooled peak | `snp_count >= 2` |
| Effect size | `abs(ref_frac - 0.5) >= 0.2` |
| ASE significance | BH-adjusted `p < 0.05` |
| Cohort caQTL or RASQUAL joint test | Study-defined permutation or FDR over the declared test family |
| Allele-frequency filter | `0.1 < AF < 0.9` |

An absolute deviation of 0.2 is a 30:70 imbalance; 0.4 is a strong 10:90
imbalance. These are source defaults, not universal biological constants.

## Failure diagnosis

### Systematic reference excess

- **Trigger:** ASEReadCounter is run on an ordinary ATAC BAM.
- **Mechanism:** Alternative alleles introduce a mismatch to the reference and
  are less likely to align.
- **Symptom:** Null sites aggregate above a 0.5 reference fraction, often around
  0.51-0.55 in the provider guidance.
- **Action:** Discard those allele counts for inference, run WASP, and recount.

### Sparse per-SNP coverage

- **Trigger:** A heterozygous site has fewer than about 30 reads.
- **Mechanism:** The binomial interval is too wide to distinguish moderate
  effects from sampling noise.
- **Action:** Pool valid replicates, aggregate consistently oriented SNPs within
  a peak, or add cohort-level evidence. Do not call non-significance balance.

### RASQUAL missing, crashing, or returning NA

- **Trigger:** A feature has no feature SNPs, or `-l` and `-m` do not match the
  SNP records in the streamed cis window.
- **Mechanism:** RASQUAL expects the reported testing-SNP and feature-SNP counts
  to describe the VCF records used for that feature.
- **Action:** Derive both counts from the actual VCF window and skip features
  with zero feature SNPs. Do not precompute or pass an external LD matrix;
  stock RASQUAL estimates genotype/allelic correlation internally.

### Unphased input used for haplotype pooling

- **Trigger:** Heterozygous genotypes are unphased.
- **Mechanism:** Reads cannot be assigned consistently to a shared haplotype
  across loci. WASP's allele-swap mapping correction and per-site GATK counts
  do not themselves make a haplotype-direction claim.
- **Action:** Phase with SHAPEIT5, BEAGLE, or whatshap before pooling loci.
  Retain a phase-set identifier and never combine separate phase blocks.

### Low cohort caQTL replication

- **Trigger:** A MatrixEQTL peak association and individual ASE disagree.
- **Mechanism:** Cohort association includes cis and trans effects, while ASE is
  within-individual cis evidence; technical or individual outliers can also
  drive disagreement.
- **Action:** Inspect individuals, run a joint model where appropriate, and
  report the disagreement instead of forcing a single conclusion.

### Multiallelic GATK failure

- **Trigger:** ASEReadCounter receives multiallelic sites.
- **Action:** Pre-filter the VCF to biallelic SNPs with bcftools.

### Slow WASP execution

- **Trigger:** Deep coverage makes realignment dominant.
- **Action:** Parallelize by chromosome or genomic region while preserving the
  same reference, indices, and merge semantics.

## Reconciliation patterns

| Pattern | Likely interpretation | Follow-up |
|---|---|---|
| GATK ASE has systematic REF > ALT | WASP was omitted or ineffective | Re-run correction and verify the null distribution |
| RASQUAL joint p-value is much smaller than its component tests | Combined total and allelic evidence increased power | Inspect all components before attributing mechanism |
| ASE appears at a SNP outside an accessibility peak | Annotation, coordinate, or non-regulatory effect may be involved | Check assembly and functional annotation |
| Cohort caQTL does not reproduce as individual ASE | Trans effect, technical artifact, or individual heterogeneity | Inspect carriers and analyze the evidence types separately |

Predicted variant effects from chromBPNet or another sequence model can be
compared with observed allelic imbalance at the same SNPs. Define the direction,
denominator, calibration set, and null before interpreting concordance. Treat it
as cross-validation, not as a replacement for mapping-bias correction or
experimental validation; this Skill does not prescribe a universal percentage.

## Literature retained from the provider

- van de Geijn B et al. 2015, *Nature Methods* 12:1061 — WASP.
- Castel SE et al. 2015, *Genome Biology* 16:195 — ASE framework and counting.
- Kumasaka N et al. 2016, *Nature Genetics* 48:206 — RASQUAL.
- Harvey CT et al. 2015, *Bioinformatics* 31:1235 — QuASAR.
- Buchkovich ML et al. 2015, *BMC Medical Genomics* 8:43 — mapping-bias
  correction with limited or absent genotype data.
- Browning SR and Browning BL. 2007, *American Journal of Human Genetics*
  81:1084 — BEAGLE phasing.
- Patterson M et al. 2015, *Journal of Computational Biology* 22:498 —
  whatshap read-based phasing.

## Related Skills

- `atac-seq/atac-peak-calling` and `atac-seq/consensus-peakset` provide peaks.
- `atac-seq/differential-accessibility` handles cohort-level accessibility.
- `atac-seq/deep-learning-atac` predicts sequence-level variant effects.
- `atac-seq/enhancer-gene-linking` maps supported variants to target genes.
- `variant-calling/vcf-basics` and `variant-calling/joint-calling` prepare VCFs.
- `phasing-imputation/haplotype-phasing` handles phasing before WASP.
- `causal-genomics/fine-mapping` uses caQTL evidence for fine mapping.
- `population-genetics/association-testing` provides the GWAS context.

