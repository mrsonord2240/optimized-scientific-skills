# Depth, frequency, missingness, soft-filter and region filters

Commands to copy; none has a wrapper. Do genotype no-calls (`routes/genotype-filter.md`) before any of the missingness steps.

```bash
# Subset by type, region or sample
bcftools view -v snps in.vcf.gz -o snps.vcf.gz
bcftools view -m2 -M2 -v snps in.vcf.gz -o biallelic_snps.vcf.gz
bcftools view -r chr1:1000000-2000000 in.vcf.gz -o region.vcf.gz
bcftools view -s sample1,sample2 in.vcf.gz -o subset.vcf.gz       # ^sample3 drops one

# Depth: inspect the distribution first, then keep the middle (extreme depth means collapsed repeats or CNV)
bcftools query -f '%DP\n' in.vcf | sort -n | \
    awk '{a[NR]=$1} END {print "5th:", a[int(NR*0.05)], "95th:", a[int(NR*0.95)]}'
bcftools filter -i 'INFO/DP>10 && INFO/DP<200' in.vcf -o depth_filtered.vcf

# Minor allele frequency (population data)
bcftools filter -i 'INFO/AF>0.01 && INFO/AF<0.99' in.vcf -o maf_filtered.vcf

# Missingness: drop bad samples, recompute per-variant call rate, then drop sites; do not apply both in one pass
bcftools stats -s - in.vcf | grep ^PSC | cut -f3,14
bcftools view -S good_samples.txt in.vcf -o sample_filtered.vcf
bcftools filter -i 'F_MISSING<0.05' sample_filtered.vcf -o site_filtered.vcf

# Soft filters: label failures at each stage, extract PASS at the end
bcftools filter -s 'LowQual' -e 'QUAL<30' in.vcf.gz | \
    bcftools filter -s 'LowDepth' -e 'INFO/DP<10' -Oz -o marked.vcf.gz
bcftools view -f PASS marked.vcf.gz -Oz -o pass_only.vcf.gz

# Exclude artifact-prone regions
bcftools view -T ^ENCFF356LFX.bed in.vcf -o filtered.vcf      # ENCODE exclusion list
bcftools view -T ^LCR-hs38.bed.gz in.vcf -o lcr_filtered.vcf   # low-complexity regions (homopolymers, simple repeats)
bcftools isec -p benchmark_dir filtered.vcf.gz truth.vcf.gz     # benchmark against GIAB
```

- `-e` excludes, `-i` includes, `-s NAME` labels instead of dropping. Operators: `< <= > >= = == != && || !`; aggregates `MIN() MAX() AVG() SUM()`.
- Guard against missing values explicitly (`INFO/DP!="."`), for the same hom-alt reason as in `routes/hard-filter.md`.
- GIAB stratification BEDs (low complexity, segmental duplications, tandem repeats) are needed for honest benchmarking.

Done when the output exists and its record count is reported against the input.
