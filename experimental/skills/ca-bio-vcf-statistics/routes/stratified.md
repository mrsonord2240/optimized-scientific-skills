# Stratified evaluation

```bash
bash scripts/stratified_stats.sh input.vcf.gz strat_out easy=easy_regions.bed difficult=difficult_regions.bed
```

Needs a bgzipped, indexed VCF and one BED per stratum (GIAB stratification BEDs: github.com/genome-in-a-bottle/genome-stratifications). Prints SNPs, indels and Ti/Tv per stratum and keeps each `bcftools stats` output. A caller at 99% overall can be at 70% in difficult regions; a Ti/Tv drop inside difficult regions confirms excess false positives there.

| Stratum | Failure mode |
|---------|--------------|
| Homopolymer runs | Systematic indel errors (Illumina, Ion Torrent) |
| Tandem repeats, low complexity | Alignment ambiguity inflates FP and FN |
| Segmental duplications | Paralogous mapping gives false hets |
| GC >70% or <25% | Coverage-bias missingness |
| MHC, centromeric | Extreme polymorphism or repetitiveness |

Done when each stratum has its Ti/Tv and counts printed.
