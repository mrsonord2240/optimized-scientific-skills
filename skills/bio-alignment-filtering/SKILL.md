---
name: bio-alignment-filtering
category: Data Analysis
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
- Python: `pysam.AlignmentFile` iteration with attribute filters (pysam; see `references/pysam.md`)
- Tag / CIGAR / read-group filters: `samtools view -e`, `-r` (`references/expressions-and-read-groups.md`); downsampling: `-s` (`references/subsampling.md`)

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
| minimap2 long-read presets (`map-ont`, `map-hifi`) | `-q 1` | `-q 60` |
| minimap2 short-read (`-x sr`) | `-q 1` | `-q 30`; `-q 60` is too strict (on real Illumina reads `-q 60` kept 58.5% of primary alignments, unique MAPQ mostly 48-59) |
| pbmm2 (PacBio) | `-q 1` | `-q 60` |

For Phred-scaled aligners (BWA, minimap2), MAPQ Q maps to ~10^(-Q/10) probability of wrong mapping. For STAR, the only emitted values are 0/1/3/255 (sentinels, not probabilities).

**There is no universal "drop ambiguous" threshold.** `-q 1` removes only MAPQ 0, which is the multi-mapper value for BWA and minimap2 alone; it leaves Bowtie2 and HISAT2 multi-mappers (MAPQ 1) and STAR multi-mappers (MAPQ 1 and 3) in.
Where the aligner writes an `NH` tag (STAR, HISAT2), `-e '[NH]==1'` selects uniquely mapped alignments on any MAPQ scale. BWA, Bowtie2 and minimap2 write no `NH`, and the expression then silently returns nothing.

Checked on a 60 kb synthetic genome with planted exact and diverged repeats (BWA 0.7.19, Bowtie2 2.5.5, HISAT2 2.2.3, minimap2 2.31, STAR 2.7.11b): the "drop ambiguous" threshold above removed 400/400 exact-repeat reads for every aligner, whereas `-q 1` kept 399/400 (Bowtie2), 374/400 (HISAT2) and 80/400 (STAR). `-e '[NH]==1'` matched `-q 255` on STAR and removed 400/400 on HISAT2. On a real STAR RNA-seq BAM, MAPQ 255 coincided with `NH==1` for all 5768 records. pbmm2 26.2.99 (`--preset CCS`) on a synthetic reference with a planted 3 kb exact repeat: 5/5 reads drawn from a unique region got MAPQ 60, 4/4 reads drawn from inside the repeat got MAPQ 0.

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
| Somatic short-variant (Mutect2, Strelka2) | `-F 1280 -q 1` | Light pre-filter only: Mutect2 applies its own MAPQ>=20, not-secondary, not-duplicate and non-chimeric-original filters (a real Mutect2 run, GATK 4.6.2.0, removed 120 of 5642 reads at MAPQ 0 and 10 through `MappingQualityReadFilter`); supplementary reads are kept by the pre-filter because they may carry real somatic SNVs |
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

The caller rationale in this table comes from caller documentation. Only the Mutect2 and HaplotypeCaller read filters were checked by running (GATK 4.6.2.0); Strelka2, DeepVariant, clair3, Sniffles, cuteSV, Manta, GRIDSS, Delly and SvABA were not run.

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

## Filter Pitfalls

Checked on the 1000G test BAM (samtools 1.24).

- **Read-level filters orphan mates.** `-q`, `-F` and `-e` judge each record alone, so one mate can go and the other stay (89 single-record templates after `-F 3332 -q 30`). `samtools fixmate` does not remove them. When the next tool needs complete pairs, keep only the names that still have two records:
  ```bash
  samtools view -F 3332 -q 30 -o filt.bam input.bam
  samtools view filt.bam | cut -f1 | sort | uniq -d > both.txt   # valid while at most two records per name remain (no supplementary)
  samtools view -N both.txt -o paired.bam filt.bam
  ```
- **Do not add `-f 2` for SV callers.** Proper-pair discards the discordant pairs they use (109 mapped records lost here); the SV rows keep `-F 1024` only.
- **`tlen` is signed.** The two mates carry `+N` and `-N`, so `-e 'tlen>=100 && tlen<=500'` keeps one mate of each pair (4677 records instead of 9346). `abs()` is not available in `-e`; test both signs: `-e '(tlen>=100 && tlen<=500) || (tlen<=-100 && tlen>=-500)'`.

## Reference Files

| Read this | When |
|-----------|------|
| `references/subsampling.md` | Downsampling to a fraction or a target read count, pair-consistent, with seeds (`samtools view -s`) |
| `references/expressions-and-read-groups.md` | Filtering on tags, CIGAR or other fields with `samtools view -e`; selecting by read group or library (`-r`, `-R`, `-l`, `-n`) |
| `references/pysam.md` | Filtering in Python: predicate filter, region and BED recipes, hash-based subsampling; `examples/filter_bam.py` is the command-line version |

Runnable code lives in `scripts/` (`filter_by_bed.py`, `subsample_pysam.py`, `match_read_count.sh`) and `examples/filter_bam.py`; invocations are in the reference files. Paths are relative to this Skill's folder.

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
| Subsample 10% (reproducible; see `references/subsampling.md`) | `view -s 42.1` |
| Standard filter | `view -F 3332 -q 30` |

## Related Skills

- sam-bam-basics - FLAG semantics, MAPQ-by-aligner, secondary vs supplementary
- alignment-sorting - Sort before/after filtering
- alignment-indexing - Required for region filtering
- alignment-amplicon-clipping - Primer clipping for amplicon panels
- duplicate-handling - Mark duplicates before filtering
- bam-statistics - Check filter effects
