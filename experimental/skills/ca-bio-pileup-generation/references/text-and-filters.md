## Output Format

Text pileup format: 3 shared columns, then 3 columns per input BAM (6 columns for one BAM, 9 for two):
```
chr20   1400300  G    6    ,,,.,,    BBB<BB
chr20   1400301  C    6    ,,,.,,    BB<BBB
```

| Column | Description |
|--------|-------------|
| 1 | Chromosome |
| 2 | Position (1-based) |
| 3 | Reference base (`N` if `-f` is missing or lacks the contig; lower case where the FASTA is soft-masked) |
| 4 | Read depth (= number of read symbols in column 5 = number of characters in column 6) |
| 5 | Read bases |
| 6 | Base qualities (ASCII - 33; BAQ-adjusted unless `-B`) |

A covered position whose bases are all filtered prints depth 0 and `*` in columns 5 and 6 (e.g. `chr22 1952 T 0 * *`); uncovered positions are absent unless `-a`/`-aa`.

### Read Bases Encoding

| Symbol | Meaning |
|--------|---------|
| `.` | Match on forward strand |
| `,` | Match on reverse strand |
| `ACGT` | Mismatch (uppercase = forward) |
| `acgt` | Mismatch (lowercase = reverse) |
| `^Q` | Start of read (Q = MAPQ + 33 as ASCII, e.g. `^]` = MAPQ 60) |
| `$` | End of read (follows the read's last base) |
| `+2AC` | Insertion of the 2 bases AC after this base (adds no depth; lowercase = reverse) |
| `-3CGT` | Deletion of the 3 reference bases CGT after this base |
| `*` | Deleted base (counted in depth; `#` with `--reverse-del` on the reverse strand) |
| `>` / `<` | Reference skip (intron, `N` in CIGAR), forward / reverse read; counted in depth |

## Quality Filtering Options

The defaults also drop data silently; they are listed in Pileup Options and Defaults below.

### Minimum Mapping Quality
```bash
samtools mpileup -f reference.fa -q 20 input.bam
```

### Minimum Base Quality
```bash
samtools mpileup -f reference.fa -Q 20 input.bam
```

### Combined Quality Filters
```bash
samtools mpileup -f reference.fa -q 20 -Q 20 input.bam
```

### Maximum Depth (Critical Trap)
```bash
# The default -d of each tool (Pileup Options and Defaults) silently truncates targeted / mt-DNA / amplicon /
# UMI-deduped data; bcftools' is far lower. Set -d explicitly in whichever you use
samtools mpileup -f reference.fa -d 0 input.bam        # no cap
samtools mpileup -f reference.fa -d 1000000 input.bam  # explicit high cap

# WRONG -- samtools mpileup writes text; bcftools call reads VCF/BCF and stops with
# "Failed to read from standard input: unknown file type"
samtools mpileup -f ref.fa in.bam | bcftools call -mv

# RIGHT -- bcftools mpileup writes the BCF that bcftools call reads; set -d explicitly
bcftools mpileup -d 1000000 -f ref.fa in.bam | bcftools call -mv
```
pysam `pileup()` also caps at `max_depth=8000`, and `max_depth=0` is **not** unlimited (still 8000); pass a large number such as `max_depth=1000000`.

## BAQ: Base Alignment Quality (Critical Default)

When `-f ref.fa` is passed, BAQ is enabled by default. BAQ Phred-scales the probability that a base is misaligned (HMM realignment against the reference over a small window) and reduces base quality near indels. Tradeoffs: slower; suppresses FP SNVs near indels; hurts indel detection sensitivity. In a **text** pileup it can hide an indel completely: with BAQ on, the bases of reads carrying a 3 bp deletion fall below `-Q 13` next to it (depth 8 -> 4, no `-3CGT` marker); `-B` shows all 8 reads and the marker.

| Flag | Behavior |
|------|----------|
| (default with `-f`) | BAQ on (computed against the reference; an existing BQ tag is reused) |
| `-B` / `--no-BAQ` | Disable BAQ -- raw qualities |
| `-E` / `--redo-BAQ` | Force recompute, ignoring existing BQ tags (cannot be combined with `-B`) |

**BAQ ON for:** short-read germline SNV (BWA, Bowtie2, HISAT2), short-read somatic SNV.

**BAQ OFF (`-B`) for:** long-read variant calling (ONT, PacBio HiFi), SV calling, RNA-seq near splice junctions, viral / amplicon, ultra-deep ctDNA from consensus reads (consensus quality already inflated), aDNA (qualities pre-rescaled by mapDamage).

### Library-Typed Flags Cheat Sheet

| Library | Flags |
|---------|-------|
| Short-read germline WGS (BWA) | `-q 20 -Q 20 -d 0` (BAQ on default) |
| Short-read tumor WGS | `-q 1 -Q 13 -d 0 -B` (low MAPQ kept; BAQ off) |
| Amplicon viral (ARTIC) | `-aa -A -d 600000 -B -Q 20` |
| Capture / exome | `-q 20 -Q 20 -d 0` |
| Long-read ONT R10.4+ | `-q 30 -Q 0 -B -d 0`; for `bcftools mpileup` add `--max-BQ 30` (its `ont` preset value) |
| PacBio HiFi | `-q 20 -Q 0 -B -d 0` |
| RNA-seq variants | `-q 20 -Q 20 -B -d 0` |
| Forensic / aDNA | `-q 0 -Q 0 -A -d 0 -B` |
