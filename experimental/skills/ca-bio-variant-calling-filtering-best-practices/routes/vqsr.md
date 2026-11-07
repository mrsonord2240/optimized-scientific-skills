# VQSR (large human cohort)

No wrapper script: this needs GATK and the truth resources. Run the commands below, once with `-mode SNP` and once with `-mode INDEL`, then merge.

Use it only for one deep WGS or about 30 or more jointly genotyped exomes. For a single exome or panel stop and use `routes/hard-filter.md`.

```bash
# SNP recalibration
gatk VariantRecalibrator \
    -R reference.fa -V cohort.vcf.gz \
    --resource:hapmap,known=false,training=true,truth=true,prior=15.0 hapmap.vcf.gz \
    --resource:omni,known=false,training=true,truth=true,prior=12.0 omni.vcf.gz \
    --resource:1000G,known=false,training=true,truth=false,prior=10.0 1000G.vcf.gz \
    --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 dbsnp.vcf.gz \
    -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
    -mode SNP \
    -O snp.recal --tranches-file snp.tranches

# Apply the truth-sensitivity tranche (99.7 keeps 99.7% of truth-set SNPs)
gatk ApplyVQSR \
    -R reference.fa -V cohort.vcf.gz \
    -mode SNP --recal-file snp.recal --tranches-file snp.tranches \
    --truth-sensitivity-filter-level 99.7 \
    -O snp.recalibrated.vcf.gz
```

- Indel run: same two commands with `-mode INDEL` and the Mills/1000G gold-indel resource. Merge the two recalibrated callsets with `bcftools concat -a` then `sort`.
- Exomes: omit `-an DP`.
- Tranches are truth-set sensitivities, not false-discovery rates.
- A converged run is not a good run; on too few variants VQSLOD is noise.

Done when both modes ran, the merged VCF exists, and `routes/validate.md` was run on it.
