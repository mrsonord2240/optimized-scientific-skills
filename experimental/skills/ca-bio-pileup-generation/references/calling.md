## Variant Calling Pipeline (Modern: bcftools mpileup)

**Goal:** Call variants from alignment data using the pileup-based approach.

**Approach:** Use `bcftools mpileup` (not `samtools mpileup -g`) so genotype-likelihood code is co-versioned with `bcftools call`. Apply quality and depth caps explicitly; annotate FORMAT fields needed for downstream filtering.

### Modern Germline Calling
```bash
bcftools mpileup -f reference.fa -d 1000000 -q 20 -Q 20 \
    --annotate FORMAT/AD,FORMAT/DP,FORMAT/SP,INFO/AD \
    input.bam | \
  bcftools call -mv -Oz -o variants.vcf.gz
bcftools index -t variants.vcf.gz
```

### BCF Intermediate (re-run `bcftools call` without repeating mpileup)
```bash
bcftools mpileup -f reference.fa -d 1000000 -Ob -o raw.bcf input.bam
bcftools call -mv raw.bcf -o variants.vcf
```

### Multi-Sample Joint Calling
```bash
bcftools mpileup -f reference.fa --threads 4 -d 100000 -q 20 -Q 20 \
    -a FORMAT/AD,FORMAT/DP s1.bam s2.bam s3.bam | \
  bcftools call -mv --threads 4 -Oz -o joint.vcf.gz
```

Joint calling over all samples beats merging per-sample VCFs. `-d` is per file: set it above the expected depth, but `-d 1000000` with several samples makes bcftools warn "Potential memory hog" (checked with 2 samples).

### Parallel by Contig
```bash
# Merge in header order: a chr*.vcf.gz glob sorts chr1, chr10, chr11, chr2 ... and bcftools concat
# accepts that silently (exit 0, valid index). Keep contigs.txt to primary contigs (no ':' or '*' names).
samtools idxstats input.bam | cut -f1 | grep -v '^\*$' > contigs.txt
xargs -a contigs.txt -P 4 -I{} bash -c \
    'set -o pipefail; bcftools mpileup -f reference.fa -r {} -d 1000000 input.bam | bcftools call -mv -Oz -o {}.vcf.gz' \
  && sed 's/$/.vcf.gz/' contigs.txt > vcf_list.txt \
  && bcftools concat -f vcf_list.txt -Oz -o all.vcf.gz && bcftools index all.vcf.gz
```

For somatic / low-VAF, prefer Mutect2 / Strelka2 / DeepVariant -- materially better than mpileup-based callers.

### Overlap Detection Defaults

When fragment length < 2 * read_length, R1 and R2 overlap. Both `samtools mpileup` and `bcftools mpileup` enable overlap detection by default and count overlapping bases once; pass `-x` to disable (long forms: `--disable-overlap-removal` / `--ignore-overlaps-removal` in samtools 1.24, `--ignore-overlaps` in bcftools). Disabling overlap correction can inflate somatic VAFs at sites covered by overlapping pairs (especially in cfDNA / FFPE).
