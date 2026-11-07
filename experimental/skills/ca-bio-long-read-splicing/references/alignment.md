## Splice-Aware Alignment

```bash
# Annotation junctions for minimap2 (BED12); --junc-bed makes it prefer annotated junctions and rescues annotated microexons (ONT needs --junc-bonus 16 for the short ones, see Microexons)
gffread gencode.v45.annotation.gtf --bed -o annotation.bed12

# PacBio HiFi (Iso-Seq) -> minimap2 splice:hq preset
minimap2 -ax splice:hq --secondary=no --junc-bed annotation.bed12 \
    -t 16 \
    reference.fa \
    isoseq.fastq.gz | \
    samtools sort -@ 8 -o isoseq_aligned.bam
samtools index isoseq_aligned.bam

# ONT direct cDNA (PCS-114, PCB-114): reads come in both orientations, so no -uf
minimap2 -ax splice -k14 --secondary=no --junc-bed annotation.bed12 --junc-bonus 16 \
    -t 16 \
    reference.fa \
    ont_cdna.fastq.gz | \
    samtools sort -@ 8 -o ont_cdna_aligned.bam
samtools index ont_cdna_aligned.bam

# ONT direct RNA (RNA004): every read is in transcript orientation, so -uf is correct
minimap2 -ax splice -uf -k14 --secondary=no --junc-bed annotation.bed12 --junc-bonus 16 \
    -t 16 \
    reference.fa \
    ont_rna.fastq.gz | \
    samtools sort -@ 8 -o ont_rna_aligned.bam
samtools index ont_rna_aligned.bam
```

**`-uf` only for reads that are all in transcript orientation** (direct RNA, orientation-fixed FLNC, stranded preps). On mixed-orientation reads it forces half of them onto the wrong strand and creates false junctions. Measured on simulated reads (pysam intron chains vs planted truth): unstranded ONT cDNA `splice -uf -k14` 47% exact chains and 50% of reads with a false junction, versus 94% and 1.1% for `splice -k14`; unoriented HiFi with `splice:hq -uf`: 29% of reads with a false junction. On oriented reads (HiFi FLNC, direct RNA) `-uf` and no `-uf` gave identical chains, so omitting it costs nothing. On real LRGASP cDNA reads, `-uf` cut the reads on annotated junctions from 92% to 54-59% and FLAIR-corrected reads from 1536 to 972-1353.

Orientation check, on an alignment made **without** `-uf` (minimap2 tags spliced reads `ts:A:+` when the read is in transcript orientation):

```bash
samtools view -F 2308 aligned.bam | awk '{for(i=12;i<=NF;i++) if($i ~ /^ts:A:/){n++; if($i=="ts:A:+") p++}} END{printf "%.3f\n", p/n}'
# ~1.000 -> oriented reads (-uf is safe); ~0.5 -> unoriented (never -uf). Measured: oriented HiFi 1.000, unstranded ONT 0.501, real LRGASP cDNA 0.503
```

`--secondary=no` discards secondary alignments. `--junc-bed` (with or without `--junc-bonus 16`) does not block novel junctions: a novel 24-nt acceptor shift stayed correctly aligned. Without an annotation, drop `--junc-bed`.

`splice:hq` is the preset for HiFi; plain `splice -k14` is for ONT. The presets matter little on simulated reads (HiFi: same chains within 2 reads; ONT with `splice:hq`: 2.1% versus 1.1% of reads with a false junction), so the main alignment risk is `-uf`, not the preset.

### Microexons (3-27 nt)

Use annotation junctions with minimap2: **HiFi `--junc-bed` alone; ONT cDNA/direct RNA `--junc-bed --junc-bonus 16`**. A high bonus can recode reads that skip the exon into inclusion reads. Test both isoforms on a known locus before raising it: inclusion may rise, but skipping must not fall; intron-next-to-soft-clip reads must remain zero. This recipe was checked across 60 simulated 4-13 nt exons (200 inclusion + 200 skipping reads each) and independently on 12 sizes from 3-27 nt with tandem pairs. Direct-RNA results vary with error profile: at bonus 16, mean inclusion/skipping was 99.1/99.2% at ~4% simulated error and 93.2/94.9% at 7.2%; do not treat one error profile as a platform guarantee. On real LRGASP ONT cDNA, bonus 16 also changed 109/1883 (5.8%) chains by extending them through annotated introns, while keeping zero dangling junctions. See [microexon validation](references/microexons.md) for the measurements, controls, and short-read-junction route.

```bash
samtools view -F 2308 aligned.bam | awk '$6 ~ /N[0-9]+S$/ || $6 ~ /^[0-9]+S[0-9]+N/' | wc -l   # reads with an intron next to a soft clip
```

```python
# usage: python check.py aligned.bam UP_START UP_END DOWN_START DOWN_END SKIP_START SKIP_END
# intron coordinates 0-based half-open (BED12 block gaps): UP/DOWN = introns flanking the microexon, SKIP = exon1 end to exon3 start
import sys, pysam
bam, *c = sys.argv[1:]; c = list(map(int, c)); up, down, skip = tuple(c[0:2]), tuple(c[2:4]), tuple(c[4:6])
inc = exc = 0
for r in pysam.AlignmentFile(bam):
    if r.is_unmapped or r.is_secondary or r.is_supplementary: continue
    pos, introns = r.reference_start, set()
    for op, n in r.cigartuples:
        if op == 3: introns.add((pos, pos + n))
        if op in (0, 2, 3, 7, 8): pos += n
    inc += up in introns and down in introns
    exc += skip in introns
print("inclusion reads", inc, "skipping reads", exc)
```

Example output for a simulated 4-nt exon (200 inclusion + 200 skipping HiFi reads): plain `0 / 399`, default bonus `199 / 201`, bonus 20 `400 / 0` (every skipping read recoded).

An unannotated microexon needs externally supported junctions, for example a 6-column BED from STAR `SJ.out.tab`. Do not substitute uLTRA as an equivalent rescue: its 150/150 result was one annotated 10-nt sequence; on the independent 3-27 nt set it lost 3-nt exons and a 4+12 nt tandem pair and was weaker below about 12 nt, especially on ONT/direct RNA. See [microexon validation](references/microexons.md).

```bash
awk 'BEGIN{OFS="\t"} $4>0 {print $1,$2-1,$3,"sj"NR,$7,($4==1?"+":"-")}' SJ.out.tab > sr_junctions.bed
# HiFi; for ONT add --junc-bonus 16
minimap2 -ax splice:hq --secondary=no --junc-bed sr_junctions.bed -t 16 reference.fa isoseq.fastq.gz | samtools sort -o isoseq_aligned.bam
```
