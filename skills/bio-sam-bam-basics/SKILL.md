---
name: bio-sam-bam-basics
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
- Python: `pysam.AlignmentFile()` (pysam)
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
`samtools view input.bam r1 r2` runs one query per region, so a record overlapping two regions is printed twice, silently (exit 0). On the test BAM `chr22:2000-3000 chr22:2500-3500` gave 7356 rows for 5426 distinct records (full-scan truth). De-duplicate with `-M` (multi-region iterator), or use a BED file:
```bash
samtools view -M input.bam chr1:1000-2000 chr1:1500-2500
samtools view -M -L regions.bed input.bam      # BED is 0-based, half-open
samtools view --region-file regions.bed input.bam
```
All three gave 5426 on the example above. In pysam, merging overlapping intervals is not enough: a read that spans the gap between two nearby intervals is fetched by both. Use:
```python
def fetch_regions(bam, regions):
    '''Yield each read overlapping any (contig, start, end) region once; 0-based, half-open.'''
    merged = {}
    for contig, start, end in sorted(regions):
        ivs = merged.setdefault(contig, [])
        if ivs and start <= ivs[-1][1]:
            ivs[-1][1] = max(ivs[-1][1], end)
        else:
            ivs.append([start, end])
    for contig, ivs in merged.items():
        prev_end = None
        for start, end in ivs:
            for read in bam.fetch(contig, start, end):
                if prev_end is None or read.reference_start >= prev_end:  # else already yielded for the previous interval
                    yield read
            prev_end = end
```

### Count Alignments
```bash
samtools view -c input.bam            # records, including secondary/supplementary
samtools view -c -F 2304 input.bam    # primary alignments only
```
Add `-@ N` to any `samtools view` for extra compression/decompression threads.

## Format Conversion

**Goal:** Convert between SAM (text), BAM (binary), and CRAM (reference-compressed) alignment formats.

**Approach:** Use `samtools view` with format flags (`-b` for BAM, `-C` for CRAM, `-h` for SAM with header). CRAM requires a reference FASTA with `-T`, for reading as well as writing. `examples/convert_formats.sh <in> <out> [reference.fa]` does this by extension (case-insensitive), passes the reference for CRAM input too, and refuses input == output.

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
CRAM stores optional tags in its own order (NM/MD move to the end), so a round trip keeps every field and tag but not the tag order.

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

## MAPQ Is Not Portable Across Aligners

`samtools view -q 30` does different things depending on what produced the BAM. MAPQ is an aligner-specific scale, not a universal probability:

| Aligner | MAPQ scale | "Unique" sentinel | Common gotcha |
|---------|-----------|-------------------|----------------|
| BWA-MEM / BWA-MEM2 | 0-60 | 60 | `-q 30` is sensible "high confidence" |
| minimap2 (DNA / pbmm2) | 0-60 | 60 | Spec-compliant |
| HISAT2 | 0-60 | 60 | Spec-compliant |
| Bowtie2 | 0-42 end-to-end (0-44 with `--local`) | 42 (44 in `--local`) is the top score and most common: 97% of records in a checked run | `-q 60` drops everything; `-q 23` is a common "uniquely mapped" convention (not a probabilistic 99% threshold) |
| STAR | 0, 1, 3, 255 | **255 = uniquely mapped (sentinel, not a quality)** | `-q 255` for "unique only"; `-q 30` accidentally keeps unique only too |
| DRAGEN | 0 to `--mapq-max` (default 60) | varies | `-q 30` still meaningful; distribution shape differs |
| Cell Ranger / STARsolo | inherits STAR | 255 | Same trap as STAR |

MAPQ 255 means "not available" in the SAM spec; only STAR (and tools that inherit it) use it for "unique". Verify the actual scale of any unfamiliar BAM:
```bash
samtools view input.bam | awk '{print $5}' | sort -un | head
samtools view -H input.bam | grep '^@PG' | head -1   # which aligner produced this BAM
```

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

## Context-Specific Tags

Beyond the standard fields, downstream tools depend on optional tags whose presence depends on aligner and assay. Inspect with `samtools view input.bam | head -1 | tr '\t' '\n'` or pysam `read.get_tag('XX')`.

| Tag | Set by | Meaning | Required by |
|-----|--------|---------|-------------|
| NM:i | bwa, minimap2, samtools calmd | Edit distance to reference | Edit-distance filters, e.g. `samtools view -e '[NM]<=2'` |
| MD:Z | bwa, samtools calmd | Mismatch positions (text) | Not needed by `bcftools mpileup` or mapDamage (output identical with MD/NM stripped, checked); `samtools calmd` regenerates it |
| MC:Z | samtools fixmate -m (bwa mem also writes it) | Mate CIGAR | samtools markdup |
| ms:i | samtools fixmate -m | Mate score (lowercase per SAMtags); minimap2's own `ms:i` is an unrelated DP score | samtools markdup |
| RG:Z | aligner from -R | Read group ID | GATK BQSR, MarkDuplicates LB lookup |
| SA:Z | All split-read aligners | Other alignments of the read: `rname,pos,strand,CIGAR,mapQ,NM;` records (pos 1-based, each ends with `;`) | Sniffles, Manta, cuteSV, GRIDSS, Delly |
| NH:i | STAR, HISAT2 | Number of reported hits | featureCounts multimapper handling, Salmon |
| HI:i | STAR | Hit index among NH (1-based by default; `--outSAMattrIHstart 0` for 0-based) | RSEM |
| XS:A | STAR (`--outSAMstrandField intronMotif`), HISAT2 | Strand inferred from splice motif | StringTie, Cufflinks |
| ts:A | minimap2 `-ax splice` | Transcript strand from splice motif | StringTie |
| CB:Z | Cell Ranger, STARsolo | Corrected cell barcode | scRNA quantification |
| UB:Z | Cell Ranger, STARsolo | Corrected UMI | UMI-aware dedup |
| RX:Z | fgbio AnnotateBamWithUmis | Raw UMI (bulk) | fgbio GroupReadsByUmi |
| MI:Z | fgbio GroupReadsByUmi | Molecular identifier (UMI group) | CallMolecularConsensusReads, duplex calling |
| cs:Z | minimap2 --cs | Compact CIGAR-with-bases | paftools, SV tools |

A missing tag can make a tool refuse or quietly do less. `samtools markdup` on a file without fixmate's tags refuses (`no ms score tag. Please run samtools fixmate on file first.`, exit 1); check each tool's tag requirements instead of assuming.

## Provenance: @PG Chain

The `@PG` lines record every tool that touched the BAM, linked through `PP` (previous program) tags. This is the audit trail.

```bash
samtools view -H input.bam | grep '^@PG'
```
`samtools view -H` appends its own `@PG` line to the output; add `--no-PG` to see the file's chain unchanged.

A clean germline pipeline:
```
@PG ID:bwa-mem PN:bwa VN:0.7.17
@PG ID:samtools.1 PN:samtools VN:1.20 PP:bwa-mem CL:samtools sort
@PG ID:samtools.2 PN:samtools VN:1.20 PP:samtools.1 CL:samtools fixmate
@PG ID:samtools.3 PN:samtools VN:1.20 PP:samtools.2 CL:samtools markdup
```

A broken/missing chain (no PP, unknown tools, gaps) means the BAM cannot be reliably reproduced.

## CRAM Reference Resolution (Critical)

CRAM stores reads relative to a reference; without it, the file is unreadable. htslib resolves the reference in this order:

1. Command-line `-T ref.fa` / `--reference`
2. `REF_CACHE` env var (local MD5-named cache; searched *before* `REF_PATH`)
3. `REF_PATH` env var (colon-separated; each element matched by the `@SQ M5:` MD5). A remote server such as EBI ENA is consulted only if its URL is present here -- it was the built-in default through htslib 1.21, but that default was **removed in 1.22** to reduce EBI load, so modern htslib does no network lookup unless that URL is added explicitly.
4. Local file named in the `@SQ UR:` header tag (local / `file://` paths only; `http`/`ftp` URIs in `UR:` are ignored)

On HPC nodes without internet, populate a local cache once:
```bash
mkdir -p $HOME/cram_cache
seq_cache_populate.pl -root $HOME/cram_cache reference.fa
export REF_CACHE=$HOME/cram_cache/%2s/%2s/%s
export REF_PATH=$REF_CACHE   # local only; no network/ENA lookup

samtools quickcheck -v file.cram                 # header + EOF only: passes with the reference missing and with a corrupt body
samtools view -o /dev/null file.cram && echo ok  # full decode: exit 1 if the reference cannot be resolved or a slice is corrupt
```
`samtools view -c`, `flagstat` and `idxstats` never decode bases, so they succeed on a CRAM whose reference is unreachable (checked on 1.24); only a full decode (`view -o /dev/null`, or `stats`) proves it.

`samtools view -C` without `-T` does not fail when no reference can be resolved: it warns (`Enabling embed_ref=2`) and embeds the read sequences. For a self-contained CRAM that decodes with no reference, use `-T ref.fa --output-fmt-option embed_ref=1`.

CRAM can be made irreversibly lossy, but the `archive` profile is NOT how: `--output-fmt-option archive` is a *lossless* maximum-compression preset (fqzcomp quality codec, name tokenization, larger slices) that does not alter bases or qualities. Irreversible loss comes instead from explicit quality **binning** (e.g. Illumina 8-bin), which must be applied deliberately and is harmful for low-coverage / somatic / forensic / archival data. Convert against the *exact* reference the BAM was aligned to (matched by `@SQ M5:`). A wrong reference is refused, not silently accepted: each slice stores its reference MD5, so decoding fails with `MD5 checksum reference mismatch` (exit 1). Bases are silently wrong only with `--input-fmt-option ignore_md5=1`.

## pysam Python Alternative

**Goal:** Read and manipulate alignment data programmatically in Python.

**Approach:** Use `pysam.AlignmentFile` to open BAM/CRAM files, iterate over reads, and access properties like coordinates, flags, CIGAR, and tags.

### Open and Iterate
```python
import pysam

with pysam.AlignmentFile('input.bam', 'rb') as bam:
    for read in bam:
        print(f'{read.query_name}\t{read.reference_name}:{read.reference_start}')
```

### Access Header
```python
with pysam.AlignmentFile('input.bam', 'rb') as bam:
    for sq in bam.header['SQ']:
        print(f'{sq["SN"]}: {sq["LN"]} bp')
```

### Read Alignment Properties
```python
with pysam.AlignmentFile('input.bam', 'rb') as bam:
    for read in bam:
        print(f'Name: {read.query_name}')
        print(f'Flag: {read.flag}')
        print(f'Chrom: {read.reference_name}')
        print(f'Pos: {read.reference_start}')  # 0-based
        print(f'MAPQ: {read.mapping_quality}')
        print(f'CIGAR: {read.cigarstring}')
        print(f'Seq: {read.query_sequence}')
        print(f'Qual: {read.query_qualities}')
        break
```

### Check Flag Properties
```python
with pysam.AlignmentFile('input.bam', 'rb') as bam:
    for read in bam:
        if read.is_paired and read.is_proper_pair:
            if read.is_reverse:
                strand = '-'
            else:
                strand = '+'
            print(f'{read.query_name} on {strand} strand')
```

### Fetch Region
```python
with pysam.AlignmentFile('input.bam', 'rb') as bam:
    for read in bam.fetch('chr1', 1000, 2000):
        print(read.query_name)
```

Reading detects SAM/BAM/CRAM from the file, so `'r'` and `'rb'` both read any of them (`bam.is_bam` / `bam.is_cram` tell which); the mode selects the format only when writing:

| Mode | Description |
|------|-------------|
| `r` / `rb` / `rc` | Read (format auto-detected) |
| `w` | Write SAM |
| `wb` | Write BAM |
| `wc` | Write CRAM (needs `reference_filename=`) |

`bam.mapped` / `bam.unmapped` come from the BAM index and are unavailable for SAM, unindexed BAM and CRAM (0 or an error); count with a scan instead, as `examples/view_bam.py <file> [limit] [reference.fa]` does.

### Convert BAM to SAM
```python
with pysam.AlignmentFile('input.bam', 'rb') as infile:
    with pysam.AlignmentFile('output.sam', 'w', header=infile.header) as outfile:
        for read in infile:
            outfile.write(read)
```

### Convert to CRAM
```python
with pysam.AlignmentFile('input.bam', 'rb') as infile:
    with pysam.AlignmentFile('output.cram', 'wc', reference_filename='reference.fa', header=infile.header) as outfile:
        for read in infile:
            outfile.write(read)
```

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
