# Alignment fields, flags, coordinates and CIGAR


```
@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:1000
@RG	ID:sample1	SM:sample1
@PG	ID:bwa	PN:bwa	VN:0.7.17
read1	0	chr1	100	60	8M	*	0	0	ACGTACGT	FFFFFFFF	NM:i:0
```
Fields are TAB-separated (this example parses with `samtools view -b`); a SAM with no `@SQ` lines needs `-t ref.fa.fai`.

Header lines start with `@`:
- `@HD` - Header metadata (version, sort order)
- `@SQ` - Reference sequence dictionary
- `@RG` - Read group information
- `@PG` - Program used to create file

Alignment fields (tab-separated):
1. QNAME - Read name
2. FLAG - Bitwise flag
3. RNAME - Reference name
4. POS - 1-based position
5. MAPQ - Mapping quality
6. CIGAR - Alignment description
7. RNEXT - Mate reference
8. PNEXT - Mate position
9. TLEN - Template length
10. SEQ - Read sequence
11. QUAL - Base qualities
12. Optional tags (NM:i:0, MD:Z:50, etc.)



| Flag | Decimal | Meaning |
|------|---------|---------|
| 0x1 | 1 | Paired |
| 0x2 | 2 | Proper pair |
| 0x4 | 4 | Unmapped |
| 0x8 | 8 | Mate unmapped |
| 0x10 | 16 | Reverse strand |
| 0x20 | 32 | Mate reverse strand |
| 0x40 | 64 | First in pair |
| 0x80 | 128 | Second in pair |
| 0x100 | 256 | Secondary alignment |
| 0x200 | 512 | Failed QC |
| 0x400 | 1024 | PCR duplicate |
| 0x800 | 2048 | Supplementary |

### Decode Flags (Bidirectional)
```bash
# Number to mnemonics
samtools flags 147
# 0x93 147 PAIRED,PROPER_PAIR,REVERSE,READ2
samtools flags 99
# 0x63 99 PAIRED,PROPER_PAIR,MREVERSE,READ1

# Mnemonics to number
samtools flags PAIRED,PROPER_PAIR,REVERSE,READ2   # 147
```

### Secondary vs Supplementary (Different Semantics)

Two different concepts that are routinely conflated:

| Bit | Name | Meaning | Filter implication |
|-----|------|---------|--------------------|
| 0x100 (256) | Secondary | An alternative candidate alignment for the same read; not the primary location | `-F 256` is correct for SNV/indel calling on short reads |
| 0x800 (2048) | Supplementary | A piece of a chimeric/split alignment (the read is split across loci) | Carries SA:Z tag; **required** by SV callers (Manta, Sniffles, cuteSV, GRIDSS, Delly) |

`-F 2304` removes both. Strip supplementary only when downstream is small-variant calling; keep supplementary for SV calling, fusion detection, or any analysis that follows split-reads.

For `-q` thresholds, MAPQ scales differ by aligner (Bowtie2 tops out at 42, STAR uses 255 for unique): see `references/mapq-by-aligner.md`.

## 0-Based vs 1-Based Coordinates (Footgun)

| Context | Coordinate system |
|---------|-------------------|
| SAM text POS | 1-based, inclusive |
| `samtools view chr1:100-200` | 1-based, closed interval |
| `samtools faidx chr1:100-200` | 1-based, closed interval |
| BAM binary internal | 0-based, half-open |
| `pysam read.reference_start` | 0-based |
| `bam.fetch('chr1', 100, 200)` | 0-based, half-open |
| BED files | 0-based, half-open |
| VCF | 1-based |
| GFF/GTF | 1-based, inclusive |

`samtools view bam chr1:100-200` and `bam.fetch('chr1', 100, 200)` return different read sets at boundaries; the pysam equivalent of the samtools region is `bam.fetch('chr1', 99, 200)`.

## CIGAR Operations

| Op | Description |
|----|-------------|
| M | Alignment match (can be mismatch) |
| I | Insertion to reference |
| D | Deletion from reference |
| N | Skipped region (introns in RNA-seq; do NOT count as covered bases) |
| S | Soft clipping (sequence in SEQ but not aligned) |
| H | Hard clipping (sequence not in SEQ) |
| = | Sequence match (explicit) |
| X | Sequence mismatch (explicit) |
| P | Padding (rare; multiple-sequence-alignment context) |

Example: `50M2I30M` = 50 bases match, 2 base insertion, 30 bases match

| Op | Consumes query (SEQ) | Consumes reference |
|----|----------------------|--------------------|
| M, =, X | yes | yes |
| I, S | yes | no |
| D, N | no | yes |
| H, P | no | no |

Aligned reference span = sum of M/D/N/=/X (pysam `reference_length`); SEQ length = sum of M/I/S/=/X (H excluded). TLEN is the signed template length: `+` on the leftmost mate, `-` on the rightmost, `0` when unavailable.

CIGAR `M` is overloaded -- it is the union of `=` and `X`. Some aligners emit `=`/`X` directly (e.g. minimap2 with `--eqx`); `bcftools mpileup` gives identical output on `=`/`X` and `M` alignments (checked on 1.24), and `samtools calmd` rebuilds MD/NM for either. `N` operations break naive coverage calculations: a 1000 bp RNA-seq read with one 50 kb intron does not cover 50 kb. Distinguish soft-clip (`S`, bases retained) from hard-clip (`H`, bases discarded -- irreversible).
