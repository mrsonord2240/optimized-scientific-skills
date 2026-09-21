---
name: bio-alignment-filtering
description: Filter alignments by flags, mapping quality, and regions using samtools view and pysam. Use when extracting specific reads, removing low-quality alignments, or subsetting to target regions.
tool_type: cli
primary_tool: samtools
license: MIT
author: GPTomics
---

## Version Compatibility

Reference examples tested with: pysam 0.22+, samtools 1.19+ (re-checked on samtools 1.24, pysam 0.24.1)

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Alignment Filtering

**"Filter my BAM file to keep only high-quality reads"** -> Select reads by FLAG bits, mapping quality, and genomic regions using samtools view or pysam.
- CLI: `samtools view` with `-F`/`-f`/`-q`/`-L` flags (samtools)
- Python: `pysam.AlignmentFile` iteration with attribute filters (pysam)

Filter alignments by flags, quality, and regions using samtools and pysam.

## Filter Flags

| Option | Description |
|--------|-------------|
| `-f FLAG` | Include reads with ALL bits set |
| `-F FLAG` | Exclude reads with ANY bits set |
| `-G FLAG` | Exclude reads with ALL bits set |
| `-q MAPQ` | Minimum mapping quality |
| `-L BED` | Include reads overlapping regions |

## Common FLAG Values

| Flag | Hex | Meaning |
|------|-----|---------|
| 1 | 0x1 | Paired |
| 2 | 0x2 | Proper pair |
| 4 | 0x4 | Unmapped |
| 8 | 0x8 | Mate unmapped |
| 16 | 0x10 | Reverse strand |
| 32 | 0x20 | Mate reverse strand |
| 64 | 0x40 | First in pair (read1) |
| 128 | 0x80 | Second in pair (read2) |
| 256 | 0x100 | Secondary alignment |
| 512 | 0x200 | Failed QC |
| 1024 | 0x400 | Duplicate |
| 2048 | 0x800 | Supplementary |

Decode a FLAG: `samtools flags 99` -> `0x63 99 PAIRED,PROPER_PAIR,MREVERSE,READ1`.

## Filter by FLAG

### Keep Only Mapped Reads
```bash
samtools view -F 4 -o mapped.bam input.bam
```

### Keep Only Unmapped Reads
```bash
samtools view -f 4 -o unmapped.bam input.bam
```

### Keep Only Properly Paired
```bash
samtools view -f 2 -o proper.bam input.bam
```

### Remove Duplicates

`-F 1024` removes only reads whose duplicate flag is **already set**. On a BAM that was never duplicate-marked it silently removes nothing (checked: 500 of 500 reads kept on a BAM with 100 unmarked duplicates). Check first, and mark if needed (see duplicate-handling):
```bash
if [ "$(samtools view -c -f 1024 input.bam)" -eq 0 ]; then
    echo "no duplicate-flagged reads; marking first" >&2
    samtools collate -O -u input.bam tmp_collate | samtools fixmate -m -u - - | \
        samtools sort -u - | samtools markdup - marked.bam
    samtools view -F 1024 -o nodup.bam marked.bam
else
    samtools view -F 1024 -o nodup.bam input.bam
fi
```

### Keep Only Primary Alignments (Remove Secondary and Supplementary)
```bash
samtools view -F 2304 -o primary.bam input.bam
```
Primary is not unique: a multi-mapped read still has one primary record (see the aligner table for uniqueness).

### Keep Read1 / Read2 Only
`-f 64` alone also returns secondary, supplementary and unmapped records of read1; add `-F 2308` for mapped primary read1.
```bash
samtools view -f 64 -F 2308 -o read1.bam input.bam
samtools view -f 128 -F 2308 -o read2.bam input.bam
```

### Forward / Reverse Strand Only
`-F 16` alone also returns unmapped reads (they carry no strand); add the unmapped bit.
```bash
samtools view -F 20 -o forward.bam input.bam      # 20 = 4 + 16
samtools view -f 16 -F 4 -o reverse.bam input.bam
```

## Filter by Mapping Quality

### Minimum MAPQ
```bash
samtools view -q 30 -o highqual.bam input.bam
```

### MAPQ and Mapped
```bash
samtools view -F 4 -q 30 -o filtered.bam input.bam
```

### Aligner-Aware MAPQ Thresholds

MAPQ scales differ by aligner; the same `-q 30` filter does different things. See sam-bam-basics for the full MAPQ-by-aligner table. Filtering recommendations:

| Aligner | "Drop ambiguous" | "High confidence" |
|---------|------------------|-------------------|
| BWA-MEM / BWA-MEM2 | `-q 1` | `-q 30` (or `-q 60` for unique only) |
| Bowtie2 | `-q 2` (multi-mappers get MAPQ 0 **or 1**) | `-q 23` (Bowtie2 MAPQ maxes at 42 end-to-end; 23 is a community "uniquely mapped" convention derived from analysis of Bowtie2's MAPQ scoring, not stated in the official manual) |
| **STAR** | `-q 255` | `-q 255` (STAR emits only MAPQ 0/1/3/255: 255 = unique, 3 = 2 loci, 1 = 3-4 loci, 0 = >4 loci; `-q 4` through `-q 255` are all equivalent, so `-q 60` is too) |
| HISAT2 | `-q 2` (multi-mappers get MAPQ 0 **or 1**) | `-q 60` |
| minimap2 (DNA, long-read) | `-q 1` | `-q 60` |
| pbmm2 (PacBio) | `-q 1` | `-q 60` |

For Phred-scaled aligners (BWA, minimap2), MAPQ Q maps to ~10^(-Q/10) probability of wrong mapping. For STAR, the only emitted values are 0/1/3/255 (sentinels, not probabilities).

**There is no universal "drop ambiguous" threshold.** `-q 1` removes only MAPQ 0, which is the multi-mapper value for BWA and minimap2 alone; it leaves Bowtie2 and HISAT2 multi-mappers (MAPQ 1) and STAR multi-mappers (MAPQ 1 and 3) in.
Where the aligner writes an `NH` tag (STAR, HISAT2), `-e '[NH]==1'` selects uniquely mapped alignments on any MAPQ scale. BWA, Bowtie2 and minimap2 write no `NH`, and the expression then silently returns nothing.

Checked on a 60 kb synthetic genome with planted exact and diverged repeats (BWA 0.7.19, Bowtie2 2.5.5, HISAT2 2.2.3, minimap2 2.31, STAR 2.7.11b): the "drop ambiguous" threshold above removed 400/400 exact-repeat reads for every aligner, whereas `-q 1` kept 399/400 (Bowtie2), 374/400 (HISAT2) and 80/400 (STAR). `-e '[NH]==1'` matched `-q 255` on STAR and removed 400/400 on HISAT2. On a real STAR RNA-seq BAM, MAPQ 255 coincided with `NH==1` for all 5768 records. The pbmm2 row was not run.

## Filter by Region

Region queries need an index (`samtools index input.bam`); a contig name absent from the header prints a warning and returns zero reads with exit code 0, so check names with `samtools idxstats`. Regions are 1-based and inclusive.

### Single Region
```bash
samtools view -o region.bam input.bam chr1:1000000-2000000
```

### Multiple Regions
```bash
samtools view -o regions.bam input.bam chr1:1000-2000 chr2:3000-4000
```

### Regions from BED File
```bash
samtools view -L targets.bed -o targets.bam input.bam
```
Each overlapping read is written once, even when it overlaps several BED rows (BED is 0-based half-open; header lines starting `track`/`browser`/`#` are skipped). Add `-P` to also retrieve the mates of reads in the region when the mates lie outside it (needs a region or `-L`, and an index).

### Combine Region and Quality
```bash
samtools view -q 30 -L targets.bed -o filtered.bam input.bam
```

## Combined Filters

### Standard Quality Filter

**Goal:** Produce a clean BAM containing only primary, mapped, non-duplicate reads with high mapping confidence.

**Approach:** Combine FLAG exclusion (-F for unmapped + secondary + duplicate + supplementary) with a MAPQ threshold. This is the one "standard" filter; duplicates must already be marked (see Remove Duplicates).

**Reference (samtools 1.19+):**
```bash
samtools view -F 3332 -q 30 -o filtered.bam input.bam
```

### Variant Calling Prep -- Assay-Aware

**Goal:** Choose a filter that matches what the downstream caller expects. Stripping supplementary alignments breaks SV callers; requiring proper-pair drops valid spliced RNA-seq reads.

| Assay / caller | Recommended filter | Why |
|----------------|-------------------|-----|
| Germline WGS short-variant (HaplotypeCaller, DeepVariant) | `-f 2 -F 3328 -q 20` | Primary, no dup, proper pair, MAPQ>=20 |
| Somatic short-variant (Mutect2, Strelka2) | `-F 1280 -q 1` | Drop only MAPQ=0; somatic callers handle low MAPQ; supplementary (chimeric) reads are kept because they may carry real somatic SNVs |
| Long-read short-variant (clair3, DeepVariant ONT) | `-F 3328 -q 5` | Long-read MAPQ scale is lower |
| Long-read SV (Sniffles, cuteSV) | `-F 1024` only | **Keep supplementary** -- SA tag is the SV signal |
| Short-read SV (Manta, GRIDSS, Delly, SvABA) | `-F 1024` only | Same -- supplementary required |
| ChIP-seq peak calling | `-F 1804 -q 30` | Drop dup + secondary + unmapped + mate-unmapped + QC-fail (supplementary reads are kept -- 1804 has no 2048 bit) |
| ATAC-seq | `-F 1804 -q 30 -f 2` | Same plus proper pair |
| Coverage analysis | `-F 1284 -q 1` | Mapped, no secondary, no dups; supplementary retained |
| RNA-seq quantification (STAR) | `-q 255` | Unique only (STAR sentinel) |
| RNA-seq quantification (HISAT2) | `-F 256 -q 60` | Different aligner semantics |
| RNA-seq variant (after `SplitNCigarReads`) | `-F 3328 -q 20` | Standard germline after split-N-trim |
| Panel / amplicon | After `samtools ampliconclip`; `-F 1024 -q 20` | Primer overlap makes proper-pair unreliable |
| ctDNA / cfDNA (UMI) | After fgbio consensus; do not pre-filter raw | |

**Reference (samtools 1.19+):**
```bash
# Short-variant germline
samtools view -f 2 -F 3328 -q 20 -o clean.bam input.bam

# SV calling: KEEP supplementary
samtools view -F 1024 -o sv_input.bam input.bam   # NOT -F 2304, 2308, 3328 or 3332

# ChIP-seq / ATAC-seq common filter
samtools view -F 1804 -q 30 -o filtered.bam input.bam
```

**Cost of getting this wrong:** filtering out supplementary reads (`-F 2304`, `2308`, `3328` or `3332`) before SV calling produces zero SV calls -- a single-flag mistake that silently invalidates the analysis.

Flag breakdowns:
- 1280 = 256 + 1024 (secondary + duplicate)
- 1284 = 4 + 256 + 1024 (unmapped + secondary + duplicate)
- 1804 = 4 + 8 + 256 + 512 + 1024 (unmapped + mate-unmapped + secondary + QC-fail + duplicate)
- 2304 = 256 + 2048 (secondary + supplementary)
- 2308 = 4 + 256 + 2048 (unmapped + secondary + supplementary)
- 3328 = 256 + 1024 + 2048 (secondary + duplicate + supplementary)
- 3332 = 4 + 256 + 1024 + 2048 (unmapped + secondary + duplicate + supplementary)

## Subsample Reads (Deterministic, Pair-Consistent)

`samtools view -s SEED.FRAC` (same as `--subsample 0.FRAC --subsample-seed SEED`) -- integer is the hash seed; fractional is the keep fraction. The hash is on QNAME, so:
1. Mate consistency: read1 and read2 are kept or dropped together.
2. Reproducibility: same input file + same seed + same fraction returns the same reads. A bare `-s 0.1` is deterministic too (it means seed 0; two runs gave byte-identical output), but always write the seed explicitly: the seed used when none is given changed in samtools 1.24 (`--subsample 0.1` without `--subsample-seed` used seed 0 before 1.24; since 1.24 it derives the seed from a hash of the input header, `--subsample-seed auto`, and picks different reads: 1001 vs 1015 records on the test BAM).
3. **Sequential downsampling requires different seeds.** With the same seed the keep-sets are nested, so `-s 1.5` then `-s 1.25` keeps a nested 1/4 (25%) of the original, not 12.5%. Use different integer seeds for independent samples.

```bash
# 10% with seed 42 (always the same reads; pair-consistent)
samtools view -s 42.1 -b -o subset.bam input.bam

# Sequential cuts with INDEPENDENT seeds
samtools view -s 1.5 -b in.bam > half1.bam
samtools view -s 2.25 -b half1.bam > quarter.bam   # 12.5% of original

# Coverage-matching to a target read count (hash-based, lands a few % off the target)
total=$(samtools view -c -F 2304 input.bam)
target=10000000
if [ "$total" -le "$target" ]; then
    echo "only $total primary reads, fewer than the target; copying unchanged" >&2
    cp input.bam matched.bam
else
    frac=$(awk -v t=$target -v n=$total 'BEGIN{printf "%.6f", t/n}')
    samtools view -s "1.${frac#*.}" -b -o matched.bam input.bam
fi

# Tumor-normal coverage matching (pull tumor down to normal)
normal_reads=$(samtools view -c -F 2308 normal.bam)
tumor_reads=$(samtools view -c -F 2308 tumor.bam)
if [ "$tumor_reads" -gt "$normal_reads" ]; then
    frac=$(awk -v n=$normal_reads -v t=$tumor_reads 'BEGIN{printf "%.6f", n/t}')
    samtools view -s "1.${frac#*.}" -b -o tumor_matched.bam tumor.bam
fi
```

The `if` guards matter: a fraction of 1 or more spliced into `-s` (`1.084419`) silently keeps only 8% of the reads.

## Expression Filtering

`samtools view -e EXPR` (or `--expr`, since samtools 1.12) supports arbitrary expression filtering on tags, FLAG, MAPQ, RNAME, CIGAR, etc. Powerful for filtering by `NM`, `AS`, `NH`, `cs`, etc. that the FLAG-based filters cannot reach (the `sclen` keyword used below is documented from samtools 1.16):
```bash
# Reads with >=2 mismatches (NM tag)
samtools view -e '[NM] >= 2' in.bam

# Soft clip on the left, on chr1
samtools view -e 'cigar=~"^[0-9]+S" && rname=="chr1"' in.bam

# Combine with FLAG and MAPQ
samtools view -F 2308 -q 30 -e '[NM] <= 5 && [AS] >= 100' in.bam

# Drop reads with low mapped fraction (samtools-internal helpers)
samtools view -e 'sclen / qlen < 0.2' in.bam
```

Note: `![NM]` is true only if NM is missing (checked on 1.24); NULL values from missing tags propagate through arithmetic. A tag the file never carries (`[XY]`, or `[NH]` on BWA output) matches nothing and gives an empty result with exit code 0 and no warning.

## Filter by Read Group
`-r` takes the read group **ID** (the `ID:` of an `@RG` header line), not a library name, and also outputs reads that carry no RG tag:
```bash
samtools view -r SRR702039 in.bam              # single read group ID (plus untagged reads)
samtools view -e '[RG]=="SRR702039"' in.bam    # strictly that read group
samtools view -R rg_list.txt in.bam            # multiple IDs via file (one ID per line)
samtools view -l LIBRARY in.bam                # by library (@RG LB: field)
```
Samtools 1.24 adds `-n` (`--exclude-no-read-group`) to drop untagged reads when `-r`/`-R` is used.

## pysam Python Alternative

`examples/filter_bam.py` wraps the filters below as a command line (`python filter_bam.py in.bam out.bam -q 30 -d -p -P -r chr1:1000-2000`): region in samtools syntax (1-based, inclusive; contig, contig:start and commas accepted), warns when `-d` finds no duplicate flags, and indexes the output when it is coordinate-sorted.

### Filter with Function

**Goal:** Apply a multi-criteria quality filter to produce clean alignments for downstream analysis.

**Approach:** Define a predicate checking mapped status, primary alignment, duplicate flag, and MAPQ; stream reads through it. Equal to `samtools view -F 3332 -q 30` (record-for-record identical on two real BAMs).

**Reference (pysam 0.22+):**
```python
import pysam

def passes_filter(read):
    if read.is_unmapped:
        return False
    if read.is_secondary or read.is_supplementary:
        return False
    if read.is_duplicate:
        return False
    if read.mapping_quality < 30:
        return False
    return True

with pysam.AlignmentFile('input.bam', 'rb') as infile:
    with pysam.AlignmentFile('filtered.bam', 'wb', header=infile.header) as outfile:
        for read in infile:
            if passes_filter(read):
                outfile.write(read)
```

### Filter by Region

pysam `fetch` uses 0-based half-open coordinates, samtools regions are 1-based inclusive: `fetch('chr1', 999999, 2000000)` returns the reads of `chr1:1000000-2000000`.
```python
import pysam

with pysam.AlignmentFile('input.bam', 'rb') as infile:
    with pysam.AlignmentFile('region.bam', 'wb', header=infile.header) as outfile:
        for read in infile.fetch('chr1', 999999, 2000000):
            outfile.write(read)
```

### Filter from BED File

**Goal:** Extract only reads overlapping target regions defined in a BED file, each read once, in coordinate order (same records as `samtools view -L`).

**Approach:** Parse BED (skipping header lines), sort and merge the intervals, fetch each merged interval, and skip reads already written through the previous interval: such a read starts before the previous interval's end. A plain loop over BED rows writes a read once per overlapped row (802 duplicated records out of 3410 on the test BAM).

**Reference (pysam 0.22+):**
```python
import pysam

def read_bed(bed_path):
    regions = []
    with open(bed_path) as f:
        for line in f:
            if not line.strip() or line.startswith(('#', 'track', 'browser')):
                continue
            parts = line.split('\t')
            regions.append((parts[0], int(parts[1]), int(parts[2])))
    return regions

with pysam.AlignmentFile('input.bam', 'rb') as infile:
    order = {name: i for i, name in enumerate(infile.references)}
    regions = sorted((r for r in read_bed('targets.bed') if r[0] in order),
                     key=lambda r: (order[r[0]], r[1]))
    merged = []
    for chrom, start, end in regions:
        if merged and merged[-1][0] == chrom and start <= merged[-1][2]:
            merged[-1][2] = max(merged[-1][2], end)
        else:
            merged.append([chrom, start, end])

    with pysam.AlignmentFile('targets.bam', 'wb', header=infile.header) as outfile:
        prev_chrom, prev_end = None, 0
        for chrom, start, end in merged:
            if chrom != prev_chrom:
                prev_end = 0
            for read in infile.fetch(chrom, start, end):
                if read.reference_start >= prev_end:
                    outfile.write(read)
            prev_chrom, prev_end = chrom, end
```

### Subsample (Pair-Consistent)

Hash on QNAME so mates stay together (a fresh `random.random()` per read drops mates inconsistently and breaks paired-end tools). Mix the seed into a real hash: `crc32(qname) ^ seed` only flips the low bits and keeps the same reads for every seed, and `crc32(f'{seed}:{qname}')` still correlates between seeds (two seeds shared 65% of their reads). This picks different reads than `samtools view -s`.
```python
import hashlib
import pysam

fraction = 0.1
seed = 42

def keep(qname):
    digest = hashlib.blake2b(f'{seed}:{qname}'.encode(), digest_size=8).digest()
    return int.from_bytes(digest, 'big') < fraction * 2**64

with pysam.AlignmentFile('input.bam', 'rb') as infile:
    with pysam.AlignmentFile('subset.bam', 'wb', header=infile.header) as outfile:
        for read in infile:
            if keep(read.query_name):
                outfile.write(read)
```

## Output Options and Checks
```bash
samtools view -b -F 4 -o output.bam input.bam                   # BAM (a .bam -o name also gives BAM)
samtools view -C -T reference.fa -F 4 -o output.cram input.bam  # CRAM (reference required)
samtools view -h -F 4 input.bam > output.sam                    # SAM with header
samtools view -c -F 3332 -q 30 input.bam                        # count only: run before writing, compare with `-c input.bam`
samtools quickcheck filtered.bam && echo OK || echo CORRUPT
samtools flagstat filtered.bam
```

## Quick Reference

| Task | samtools command |
|------|------------------|
| Mapped only | `view -F 4` |
| Unmapped only | `view -f 4` |
| Properly paired | `view -f 2` |
| Primary alignments (multi-mappers still included) | `view -F 2304` |
| No duplicates (marked first) | `view -F 1024` |
| High MAPQ | `view -q 30` |
| Region | `view file.bam chr1:1-1000` |
| BED regions | `view -L file.bed` |
| Subsample 10% (reproducible) | `view -s 42.1` |
| Standard filter | `view -F 3332 -q 30` |

## Related Skills

- sam-bam-basics - FLAG semantics, MAPQ-by-aligner, secondary vs supplementary
- alignment-sorting - Sort before/after filtering
- alignment-indexing - Required for region filtering
- alignment-amplicon-clipping - Primer clipping for amplicon panels
- duplicate-handling - Mark duplicates before filtering
- bam-statistics - Check filter effects
