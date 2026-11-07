# Hard filters in Python (cyvcf2)

```bash
python scripts/cyvcf2_filter.py in.vcf.gz out.vcf min_dp=10
```

Applies the same thresholds as `routes/hard-filter.md`, SNPs and indels separately (SNP limits for SNPs, indel limits for indels). A missing annotation passes, so hom-alt sites with undefined RankSums survive. Records that are neither SNP nor indel are counted, not written, and reported on stderr.

- Unlike `hard_filter.sh`, this script has a site DP floor (`min_dp=10`, a missing DP passes). `min_dp=0` removes it and makes the result match the bcftools script.
- `variant.gt_types` codes: 0 HOM_REF, 1 HET, 2 UNKNOWN (missing), 3 HOM_ALT.
- cyvcf2 has no Windows wheel; run on Linux, macOS or WSL.

Done when the script printed its counts dictionary and the output VCF exists.
