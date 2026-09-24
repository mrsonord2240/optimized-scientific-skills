---
name: bio-sam-bam-basics
category: Data Analysis
description: View, convert, and understand SAM/BAM/CRAM alignment files using samtools and pysam. Use when inspecting alignments, converting between formats, or understanding alignment file structure.
tool_type: cli
primary_tool: samtools
license: MIT
author: GPTomics
---

## Version Compatibility

Reference examples tested with: pysam 0.22+, samtools 1.19+; behaviour re-checked on samtools 1.24, pysam 0.24.1, bcftools 1.24

Install: `conda install -c bioconda samtools pysam` (or `pip install pysam`; pysam has no Windows wheel, use WSL/Linux).

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# SAM/BAM/CRAM Basics

**"Read a BAM file"** -> Open a binary alignment file and iterate over aligned reads with their mapping coordinates, flags, and quality scores.
- Python: `pysam.AlignmentFile()` (pysam; see `references/pysam.md`)
- CLI: `samtools view` (samtools)
- R: `scanBam()` (Rsamtools)

View and convert alignment files using samtools and pysam.

## Format Overview

| Format | Description | Use Case |
|--------|-------------|----------|
| SAM | Text format, human-readable | Debugging, small files |
| BAM | Binary compressed SAM | Standard storage format |
| CRAM | Reference-based compression | Long-term archival, smaller than BAM |

## SAM Format Structure

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

## samtools view

### View BAM as SAM
```bash
samtools view input.bam | head
```

### View with Header
```bash
samtools view -h input.bam | head -100   # keep -h when piping to another samtools/BAM writer, or the header is lost
```

### View Header Only
```bash
samtools view -H input.bam
```

### View Specific Region
```bash
samtools view input.bam chr1:1000-2000
samtools view input.bam chr1             # whole chromosome
```
Region queries need a coordinate-sorted, indexed file (`@HD SO:coordinate`, `samtools index`); without an index: `Could not retrieve index file`. Contig names must match `@SQ SN:` exactly (`chr1` vs `1`): an unknown name prints `[main_samview] region "22:2000-3000" specifies an invalid region or unknown reference. Continue anyway.` and returns 0 records with **exit code 0**. List the names with `samtools idxstats input.bam | cut -f1`.

### Multiple Regions (overlaps print records twice)
`samtools view input.bam r1 r2` runs one query per region, so a record overlapping two regions is printed twice, silently (exit 0). On the test BAM `chr22:2000-3000 chr22:2500-3500` gave 7356 rows for 5426 distinct records (full-scan truth). De-duplicate with `-M` (multi-region iterator; needs an index, like any region query), or use a BED file:
```bash
samtools view -M input.bam chr1:1000-2000 chr1:1500-2500
samtools view -M -L regions.bed input.bam      # BED is 0-based, half-open
samtools view --region-file regions.bed input.bam
```
All three gave 5426 on the example above. In pysam, merging overlapping intervals is not enough: a read that spans the gap between two nearby intervals is fetched by both. Use `scripts/fetch_regions.py` (each read once; 0-based half-open regions):
```bash
python scripts/fetch_regions.py input.bam regions.bed          # prints SAM lines; BED as above
```
```python
from fetch_regions import fetch_regions    # scripts/ on sys.path; regions = [(contig, start, end), ...]
for read in fetch_regions(bam, regions): ...
```

### Count Alignments
```bash
samtools view -c input.bam            # records, including secondary/supplementary
samtools view -c -F 2304 input.bam    # primary alignments only
```
Add `-@ N` to any `samtools view` for extra compression/decompression threads.

## Format Conversion

**Goal:** Convert between SAM (text), BAM (binary), and CRAM (reference-compressed) alignment formats.

**Approach:** Use `samtools view` with format flags (`-b` for BAM, `-C` for CRAM, `-h` for SAM with header). CRAM requires a reference FASTA with `-T`, for reading as well as writing (resolution order, offline cache, checks: `references/cram-reference.md`). `examples/convert_formats.sh <in> <out> [reference.fa]` does this by extension (case-insensitive), passes the reference for CRAM input too, and refuses input == output.

### BAM to SAM
```bash
samtools view -h -o output.sam input.bam
```

### SAM to BAM
```bash
samtools view -b -o output.bam input.sam
```

### BAM to CRAM
```bash
samtools view -C -T reference.fa -o output.cram input.bam
```

### CRAM to BAM
```bash
samtools view -b -T reference.fa -o output.bam input.cram
```

### Pipe Conversion
```bash
samtools view -b input.sam > output.bam
```

### CRAM round trip is not byte-lossless
BAM -> CRAM -> BAM keeps read names, positions, SEQ and QUAL, but not the rest verbatim (checked on samtools 1.24): `=`/`X` CIGAR ops are rewritten to `M`; NM/MD are regenerated when the reference is available (MD:Z appears on reads that had none); tags are reordered (NM/MD move to the end); an unmapped read's MAPQ becomes 0. Use BAM when the exact CIGAR ops or tag set matter.

## Common Flags

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

## Reference Files

Read the file when the task needs it; SKILL.md above covers viewing, conversion, FLAG, coordinates and CIGAR.

| File | Read when |
|------|-----------|
| `references/mapq-by-aligner.md` | Choosing or interpreting a `samtools view -q` threshold, or a BAM of unknown origin (MAPQ scale per aligner) |
| `references/tags-and-provenance.md` | A tool needs an optional tag (NM, MD, MC, ms, SA, NH, HI, CB, UB, RX, MI...), or the `@PG` chain must be audited |
| `references/cram-reference.md` | Any CRAM read or write: reference resolution order, offline `REF_CACHE`, proving a CRAM readable, embedded references, lossy vs lossless, wrong reference |
| `references/pysam.md` | Reading or writing alignments from Python: `AlignmentFile` modes, properties, fetch, format conversion, `bam.mapped` limits |

## Quick Reference

| Task | samtools | pysam |
|------|----------|-------|
| View BAM | `samtools view file.bam` | `AlignmentFile('file.bam', 'rb')` |
| View header | `samtools view -H file.bam` | `bam.header` |
| Count records | `samtools view -c file.bam` | `sum(1 for _ in bam)` (`bam.count(until_eof=True)` for unindexed) |
| Get region | `samtools view file.bam chr1:1-1000` | `bam.fetch('chr1', 0, 1000)` |

## Related Skills

- alignment-indexing - Create indices for random access (required for fetch/region queries)
- alignment-sorting - Sort alignments by coordinate or name
- alignment-filtering - Filter alignments by flags, quality, regions
- alignment-validation - Sequence dictionary cross-validation (M5 checksums)
- bam-statistics - Generate statistics from alignment files
- reference-operations - REF_PATH/REF_CACHE setup for CRAM
- sequence-io/read-sequences - Parse FASTA/FASTQ input files
