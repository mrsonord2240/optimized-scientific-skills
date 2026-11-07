# ASB methods, prerequisites, and command surfaces

Use this reference after choosing a method in `SKILL.md`. Verify each installed tool's current help before executing these provider-supplied command patterns.

## Dependency clues

### WASP path

- WASP 0.3.4+ and Python 3
- `pysam` and PyTables (`tables`)
- the WASP SNP-table preparation utilities (`snp2h5` or the text-SNP extraction path)
- bowtie2 or another supported aligner
- samtools and bcftools
- a sample-specific heterozygous VCF

Provider installation sketch:

```bash
git clone https://github.com/bmvdgeijn/WASP.git
conda install -c bioconda bowtie2 samtools bcftools gatk4
```

### BaalChIP path

- R 4.5 and Bioconductor 3.22
- BaalChIP 1.36.0 (the tested API contract)
- GenomicRanges and rtracklayer
- an official-format sample-sheet TSV, group-named heterozygous-SNP table, peaks, blacklist, and either per-variant `RAF` or group-named gDNA BAMs

```r
BiocManager::install(version = '3.22')
BiocManager::install(c('BaalChIP', 'AllelicImbalance',
                       'GenomicRanges', 'rtracklayer'))
stopifnot(as.character(packageVersion('BaalChIP')) == '1.36.0')
```

### RASQUAL path

- RASQUAL 1.1+ built from source
- tabix/bgzip
- phased genotype VCF
- RASQUAL binary phenotype/covariate inputs produced by its `txt2bin` utilities

```bash
git clone https://github.com/dg13/rasqual.git
make -C rasqual
```

### AlleleSeq path

- [`trgaleev/AlleleSeq2`](https://github.com/trgaleev/AlleleSeq2) commit `cfe8acf` plus its legacy Python 2, STAR, Picard, GNU Make, Java, reference/liftOver resources, and the official `vcf2diploid.jar`
- a phased sample VCF and matching reference FASTA

## Method taxonomy

| Method | Approach | Strength | Failure boundary |
|---|---|---|---|
| WASP | Swap alleles in variant-overlapping reads, remap, and discard inconsistent mappings | Aligner-agnostic first-stage mapping-bias control | Reduces usable depth; misses variants absent from its SNP inputs |
| RASQUAL | Joint genotype-phenotype association with per-feature `phi` bias parameter | Combines cis-QTL and ASB evidence in cohorts | Computationally intensive; sparse features and CN imbalance challenge the model |
| BaalChIP | Bayesian beta-binomial model with explicit RM/RAF corrections | Replicate-aware inference; supports measured RAF or matching gDNA | Slower; inputs and correction behavior must match the installed release |
| AlleleSeq | Align to personalized diploid genomes | Avoids a single-reference alignment target | Requires accurate phasing and per-sample genome construction |
| MBASED | Meta-analysis-based allele-specific expression | Gene-level aggregation | RNA-oriented and less precise for narrow TF peaks |
| AllelicImbalance | Bioconductor count and test workflow | Convenient R workflow | Requires variants and BAM; not a substitute for mapping-bias correction |
| deepSEA/chromBPNet | Sequence-model variant-effect prediction | Does not require sample chromatin reads | Prediction, not an ASB measurement |

For cancer data, BaalChIP 1.36.0 supports relative-allele-frequency correction through a `RAF` column in the group het table or through group-named gDNA BAMs supplied as `CorrectWithgDNA`. A detached copy-number BED is not a supported constructor or `getASB` argument, and ordinary population `AF` is not `RAF`. Record observed WASP attrition and the correction source actually used.

## WASP mapping-bias filter

The provider workflow aligns reads, identifies reads overlapping heterozygous SNPs, swaps the alleles, remaps those reads, and retains only consistently mapped fragments.

```bash
# 1. Initial alignment
bowtie2 -x hg38 -1 R1.fq -2 R2.fq -S step1.sam
samtools view -bS step1.sam | samtools sort -o step1.bam
samtools index step1.bam

# 2. Identify SNP-overlapping reads and generate allele-swapped reads
python /path/to/WASP/mapping/find_intersecting_snps.py \
    --is_paired_end \
    --is_sorted \
    --output_dir wasp_out/ \
    --snp_tab snps_tab.h5 \
    --snp_index snps_index.h5 \
    --haplotype haplotypes.h5 \
    --samples sample_list.txt \
    step1.bam

# 3. Remap both swapped mates with exactly the initial alignment settings.
# WASP v0.3.4 emits separate fq1/fq2 files for paired input.
bowtie2 -x hg38 \
  -1 wasp_out/step1.remap.fq1.gz \
  -2 wasp_out/step1.remap.fq2.gz | \
  samtools view -b -o step2.remap.bam

# 4. Retain reads that remap consistently
python /path/to/WASP/mapping/filter_remapped_reads.py \
    wasp_out/step1.to.remap.bam step2.remap.bam step2.remap.keep.bam

# 5. Final bias-filtered BAM
samtools merge -f step1.merged.bam \
  wasp_out/step1.keep.bam step2.remap.keep.bam
samtools sort -o step1.wasp.bam step1.merged.bam
samtools index step1.wasp.bam
```

Assert that the two remap FASTQs contain the same positive fragment count, the final BAM is wholly paired, and `final = direct_keep + remap_keep <= input`. The exact names above are tested against WASP v0.3.4; inspect the installed release before adapting them.

## BaalChIP workflow

The shipped runner targets the official BaalChIP 1.36.0 contract from Bioconductor 3.22. Its sample sheet is a TSV filename with exactly the provider-required fields (additional fields are not needed):

```text
group_name  target  replicate_number  bam_name  bed_name
TUMOR       FOXA1   1                 rep1.wasp.bam  rep1.narrowPeak
TUMOR       FOXA1   2                 rep2.wasp.bam  rep2.narrowPeak
```

The het table must contain unique `ID`, `CHROM`, one-based `POS`, single-base uppercase `REF`/`ALT`, and optionally measured `RAF`. `AF` is rejected rather than relabeled. Choose exactly one correction:

- `--correction raf`: every retained variant has a finite measured reference allele frequency in `[0,1]`;
- `--correction gdna --gdna-bam matched-normal.bam`: BaalChIP derives RAF from matching genomic-DNA reads. The runner removes a supplied RAF column because official BaalChIP gives it priority over gDNA.

```bash
Rscript scripts/baalchip_workflow.R \
  --samples samples.tsv --hets tumor-hets.tsv --group TUMOR \
  --blacklist hg38-blacklist.bed --imprinted imprinted-hg38.bed \
  --sex female --assembly GRCh38 --samples-assembly GRCh38 \
  --hets-assembly GRCh38 --blacklist-assembly GRCh38 \
  --imprinted-assembly GRCh38 --correction gdna \
  --gdna-bam matched-normal.bam --samtools /opt/conda/bin/samtools \
  --out results/TUMOR
```

All assembly declarations and contig styles must agree. Required files, BAM indexes, BAM-header contigs, unique replicate numbers, and group names are checked before package loading. The runner records imprinted/blacklist/chrX exclusions and reason codes, refuses an existing output path, stages output atomically, records the selected correction source, and emits a structured `no_variants_after_filters` or `no_calls` state instead of indexing an absent report.

`BaalChIP.report` returns a named list, one data frame per group. The runner selects `report[[group]]` and verifies the documented 15 report fields before export; it never treats the list itself as a data frame. Use `--preflight-only` to validate contracts without claiming model execution.

## RASQUAL joint cis-QTL and ASB

The genotype VCF is piped from tabix; there is no `--vcf` flag. Options are single-dash. Inputs `-y`, `-k`, and `-x` are native-double binary files built by the official source conversion script, not renamed text files.

```bash
python scripts/rasqual_cohort.py prepare \
  --rasqual-source /opt/rasqual --r /opt/R/bin/R \
  --y Y.txt --k K.txt --x X.txt --out cohort-binaries
```

Create a tab-separated manifest with exactly these columns and exactly one row for every Y/K feature index:

```text
feature_id  region  feature_index  n_samples  n_test_snps  n_feature_snps  exon_starts  exon_ends  lead_only
C11orf21    11:2315000-2340000  1  24  378  62  2316875,2320655  2319151,2320937  true
```

Then run the cohort family:

```bash
python scripts/rasqual_cohort.py run \
  --manifest features.tsv --feature-count 2 --covariates 4 \
  --y cohort-binaries/Y.bin --k cohort-binaries/K.bin --x cohort-binaries/X.bin \
  --vcf genotypes.vcf.gz --tabix /usr/bin/tabix --rasqual /opt/rasqual/src/rasqual \
  --output rasqual-cohort.tsv
```

The driver validates binary sizes from the declared dimensions; requires feature indices `1..N` exactly once; invokes each feature once; checks the documented 25-column schema, feature ID, finite `phi`, and convergence status `0`; derives a one-degree-of-freedom p-value from the likelihood-ratio chi-square; and adds Benjamini-Hochberg q-values over every result row emitted by this run. The original RASQUAL within-feature column 10 is retained separately. Define the family before running; do not merge selectively retained rows afterward.

## AlleleSeq external-only route

This skill does not advertise a locally executable AlleleSeq pipeline. The tested environment could dry-run the `PIPELINE.mk` from [`trgaleev/AlleleSeq2`](https://github.com/trgaleev/AlleleSeq2) at commit `cfe8acf`, but it lacked Python 2, STAR, Picard, required reference/liftOver resources, and the official `vcf2diploid.jar`; the canonical AlleleSeq host was unavailable during the bounded tooling run. No third-party JAR may be substituted.

To restore this route, provide the official artifact and a separately isolated legacy environment containing every prerequisite named above. Populate the provider Make variables for `READS_R1`, `READS_R2`, sample `PREFIX`, reference resources, and output root, then inspect a complete `make -n` plan before execution. The provider plan performs its own STAR-based alignment; do not prepend disconnected maternal/paternal Bowtie2 SAMs and imply the Make target consumes them. Until that exact plan and toolchain are validated, route AlleleSeq requests to an external maintained installation and classify execution as unavailable rather than partial success.

## Interval filters

There is no single canonical hosted hg38 BED for imprinted genes. Derive and version an interval file from a documented catalog, map it to the analysis assembly, and record the source.

```bash
bedtools intersect -v \
  -a hetSNPs.bed \
  -b imprinted_loci_hg38.bed \
  > hetSNPs.non_imprinted.bed

awk '$1 != "chrX"' hetSNPs.bed > hetSNPs.autosomal.bed
```

For cancer samples, use ASCAT, Sequenza, FACETS, or another justified caller to identify altered segments for exclusion or stratification. Do not pass a segment BED to BaalChIP 1.36.0 as though it were consumed: the tested correction interfaces are per-variant `RAF` and `CorrectWithgDNA`.

## Claim-level provenance and exact bindings

| Claim or executable surface | Canonical source | Tested binding |
|---|---|---|
| WASP swaps variant-overlapping alleles and filters inconsistently remapped reads; paired runs emit mate-specific remap FASTQs | [van de Geijn et al., Nature Methods](https://doi.org/10.1038/nmeth.3582); [official WASP repository](https://github.com/bmvdgeijn/WASP) | tag v0.3.4, commit `f980683` |
| RASQUAL binary construction, single-dash CLI, 25 output columns, `phi`, and convergence column | [Kumasaka et al., Nature Genetics](https://doi.org/10.1038/ng.3366); [official RASQUAL repository](https://github.com/dg13/rasqual) | commit `5aa553c` |
| BaalChIP sample/het schemas, `CorrectWithgDNA`, `RAFcorrection`, and named-list report | [de Santiago et al., Genome Biology](https://doi.org/10.1186/s13059-017-1165-7); [Bioconductor 3.22 BaalChIP](https://bioconductor.org/packages/3.22/bioc/html/BaalChIP.html) | BaalChIP 1.36.0, Bioconductor 3.22, R 4.5 |
| AlleleSeq personalized-genome method and external-only legacy implementation boundary | Canonical method: [Rozowsky et al., Molecular Systems Biology](https://doi.org/10.1038/msb.2011.54); tested implementation: [`trgaleev/AlleleSeq2`](https://github.com/trgaleev/AlleleSeq2) | `trgaleev/AlleleSeq2` commit `cfe8acf`; local execution unavailable |

The version bindings above are tested contracts, not open-ended minimum-version promises. Reinspect help, schemas, and output fields before using a different release.

## Suggested user requests

- Apply WASP to a ChIP-seq BAM using the sample's heterozygous VCF before testing ASB.
- Run BaalChIP on cancer ChIP-seq while accounting for copy-number imbalance.
- Use RASQUAL for joint cis-QTL and ASB analysis in a phased population cohort.
- Route an AlleleSeq personalized-genome request to a restored official legacy environment and verify its complete Make plan before execution.
- Filter imprinted loci and handle chrX before reporting ASB.
- Compare chromBPNet predictions with measured ASB while retaining coverage and sample-context caveats; the provider prompt uses `|log2_fc| > 1` as an example prediction threshold, not a universal ASB criterion.
