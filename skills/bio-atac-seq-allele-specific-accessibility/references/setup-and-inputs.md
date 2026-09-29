# Setup and inputs

Use this reference when preparing an environment or validating files before
running allele-specific accessibility analysis.

## Installation clues

```bash
# WASP for mapping-bias correction
git clone https://github.com/bmvdgeijn/WASP

# GATK 4 for allele counting
conda install -c bioconda gatk4

# RASQUAL for joint caQTL modeling
git clone https://github.com/natsuhiko/rasqual

# Phasing and supporting command-line tools
conda install -c bioconda shapeit5 beagle whatshap
conda install -c bioconda samtools bcftools vcftools plink2 bowtie2 bwa-mem2
```

These optional installs belong to separately specified external routes; the
current Skill does not provide an executable contract for them:

```r
install.packages('MatrixEQTL')
remotes::install_github('piquelab/QuASAR')
```

QuASAR is installed from GitHub rather than CRAN or Bioconductor. WASP and
RASQUAL are source checkouts in the provider examples. Pin the checkout commit
in a reproducible analysis.

The within-peak helper additionally requires Python 3, pandas, scipy,
pybedtools, and a working bedtools installation.

## Required inputs

- A deduplicated, mapping-quality-filtered ATAC-seq BAM.
- A genotype VCF containing the single BAM read-group sample. The runnable
  pipeline selects that sample before filtering to biallelic heterozygous SNPs.
- The matching reference-genome FASTA and aligner index.
- For within-peak pooling, phased haplotypes with a non-missing phase-set (`PS`)
  value from SHAPEIT5, BEAGLE, whatshap, or an equivalent method.
- A consensus peak BED on the same genome assembly.

For the shipped RASQUAL cohort workflow, also prepare its binary count/offset
matrices, six-column feature table, sample count, and indexed cohort VCF. A
different cohort model requires its own explicit matrix, covariate, model, QC,
and output contract.

## Preflight checks

1. Confirm that BAM, VCF, FASTA, aligner index, and peaks use the same genome
   assembly and chromosome naming.
2. Confirm that the BAM read groups and the VCF sample column identify the same
   individual.
3. Ensure the BAM is coordinate-sorted and indexed.
4. Ensure the FASTA is indexed with `samtools faidx` and the VCF is bgzip
   compressed and tabix indexed.
5. Confirm that pooled sites have phased genotypes and non-missing phase-set
   identifiers. Per-site counting can use unphased heterozygotes, but unphased
   sites do not support a shared haplotype direction and must not be pooled.
6. Remove multiallelic sites before GATK ASEReadCounter.
7. Check that the BAM basename used by the example pipeline matches the sample
   identifier encoded in WASP's haplotype HDF5 inputs. A mismatch can yield no
   intersecting heterozygous sites.
8. Inspect `<tool> --help` for every installed CLI because WASP and RASQUAL
   interfaces vary across source checkouts.

The runnable pipeline accepts these positional arguments:

```text
scripts/wasp_ase_pipeline.sh \
  BAM VCF GENOME_FASTA BOWTIE2_INDEX PEAKS_BED WASP_DIR OUTPUT_DIR
```

All seven arguments are required and every path is quoted by the pipeline. The
output directory must not already exist: a run is assembled in a temporary
sibling directory and renamed only after all postconditions pass. To rerun,
choose a new output path or explicitly archive/remove the old result.

