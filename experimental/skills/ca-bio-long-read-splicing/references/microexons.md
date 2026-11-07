# Microexon validation details

The primary recipe is minimap2 with annotation junctions: HiFi uses `--junc-bed`; ONT cDNA and direct RNA use `--junc-bed --junc-bonus 16`. It was evaluated with minimap2 2.31 using CIGAR intron chains, not only inclusion counts.

On 60 simulated annotated 4-13 nt microexons (three random genomes; 200 inclusion and 200 skipping reads per gene/platform), inclusion/skipping percentages were:

| Alignment | HiFi | ONT cDNA | ONT direct RNA |
|---|---|---|---|
| plain preset | 25 / 100 | 0 / 98 | 0 / 96 |
| `--junc-bed` (bonus 9) | 100 / 100 | 49 / 100 | 42 / 100 |
| bonus 16 | 100 / 100 | 99.5 / 100 | 99.1 / 100 |
| bonus 20 | 100 / 57 | 100 / 91 | 100 / 90 |

An independent 3-27 nt set (12 single exons plus two tandem pairs; 150 reads per state) confirmed the two-way control: HiFi with `--junc-bed` gave 100.0% inclusion/99.9% skipping and ONT bonus 16 gave 99.6%/99.2%. Direct RNA gave 99.1%/99.2% at about 4% simulated error, but 93.2%/94.9% at 7.2% error; four sizes fell below 90% inclusion in the latter. On real LRGASP WTC-11 ONT cDNA (1,883 reads), bonus 16 had no dangling junctions but changed 109 chains (5.8%) to fully annotated chains. Bonus 17 produced 12 dangling junctions and caused `flair correct` to fail; use 16 as a starting point, not a universal optimum.

To validate a candidate bonus, measure both the flanking-junction inclusion chain and the exon-skipping chain, then reject a setting if skipping falls as inclusion rises. Also require zero dangling junctions:

```bash
samtools view -F 2308 aligned.bam | awk '$6 ~ /N[0-9]+S$/ || $6 ~ /^[0-9]+S[0-9]+N/' | wc -l
```

`uLTRA` is not an interchangeable rescue route. Although it recovered one annotated 10-nt test exon, it lost 3-nt exons and a 4+12 nt tandem pair in the independent set, with mean inclusion 91.4% (HiFi), 76.6% (ONT), and 78.6% (direct RNA). For an exon absent from the annotation, supply orthogonal junctions instead:

```bash
awk 'BEGIN{OFS="\t"} $4>0 {print $1,$2-1,$3,"sj"NR,$7,($4==1?"+":"-")}' SJ.out.tab > sr_junctions.bed
minimap2 -ax splice:hq --secondary=no --junc-bed sr_junctions.bed -t 16 reference.fa isoseq.fastq.gz | samtools sort -o isoseq_aligned.bam
```
