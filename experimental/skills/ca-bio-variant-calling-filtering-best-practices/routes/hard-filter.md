# Germline hard filters (SNPs and indels)

```bash
ASSAY=wgs bash scripts/hard_filter.sh in.vcf.gz out_prefix
```

`ASSAY` (`wgs` or `wes`) is optional and only sets the Ti/Tv floor for the warning. Output is `out_prefix_all_filtered.vcf.gz` with an index, plus separate SNP and indel files. The script splits by type, filters each, merges with `concat -a` and `sort`, and prints record counts and Ti/Tv before and after. It removes failing records; it does not label them.

Thresholds (GATK lenient starting points; a record fails if any applies; QUAL must be at least 30 for both types, which is not part of GATK's canonical set):

| Annotation | SNP | Indel |
|------------|-----|-------|
| QD | < 2.0 | < 2.0 |
| FS | > 60 | > 200 |
| SOR | > 3.0 | > 10.0 (a community/WDL convention, not GATK canonical) |
| MQ | < 40 | not used |
| MQRankSum | < -12.5 | not used |
| ReadPosRankSum | < -8.0 | < -20.0 |

- The SNP/indel differences are deliberate: indels sit in messier alignments (higher legitimate strand bias), spurious indels cluster at read ends, and mapping-quality metrics are weaker for indels.
- The script refuses (exit 2) a DeepVariant or DRAGEN header, a VCF missing any of QD, FS, SOR, MQ, MQRankSum, ReadPosRankSum in the INFO header, and a VCF where SNP plus indel records do not equal all records (MNPs, mixed or multiallelic SNP+indel records). For the last, run `bcftools norm -m- --atomize` first. `ALLOW_DL_CALLER=1` and `ALLOW_DROP=1` override.
- If the request names its own thresholds, edit the numbers in the script's two `bcftools filter` expressions and keep the `|| INFO/X = "."` guard on every term.
- Look at the annotation histograms (known true variants against likely false positives) before trusting defaults; GATK values are population-level starting points.
- Depth alone is a poor filter: it removes real variants in collapsed duplications. Combine DP with MQ and MQRankSum.

Done when the script printed `kept:` counts and Ti/Tv, and hom-alt sites are still present in the output.
