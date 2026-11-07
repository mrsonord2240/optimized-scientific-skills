# Contamination signals

```bash
bash scripts/het_allele_balance.sh input.vcf.gz
```

One line per sample: mean allele balance of biallelic 0/1 hets from `FORMAT/AD`, expected ~0.5. Handles 0/1 and 1/0 (phased or not) only.

- A shift off 0.5 together with a raised het count and het/hom ratio, and a lower novel-fraction Ti/Tv, is a contamination signal, not the measurement.
- The measurement runs on the BAM/CRAM: VerifyBamID2 estimates the contamination fraction alpha; CHARR estimates it from a VCF or gVCF. Alpha above ~0.02-0.03 is a red flag; somatic pipelines are sensitive to 1%.
- GATK pipelines take it as `--contamination alpha`.

Done when per-sample allele balance is read, and any flagged sample is sent for VerifyBamID2 or CHARR.
