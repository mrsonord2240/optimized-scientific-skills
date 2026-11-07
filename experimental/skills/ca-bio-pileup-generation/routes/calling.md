# Pileup-based variant calling

```bash
bcftools mpileup -f reference.fa -d 1000000 -q 20 -Q 20 -a FORMAT/AD,FORMAT/DP,FORMAT/SP,INFO/AD input.bam | bcftools call -mv -Oz -o variants.vcf.gz
```

Index the FASTA with `samtools faidx reference.fa`. Set the per-file depth cap for expected coverage. Index the result with `bcftools index -t variants.vcf.gz`. Text from samtools mpileup cannot feed bcftools call. See `references/calling.md` for BCF intermediates, joint calling, parallel contigs and alternative callers, and `references/text-and-filters.md` for assay-specific filters.

Done when the requested BCF or VCF is written and its settings are recorded.
