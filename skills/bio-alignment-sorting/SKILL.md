---
name: bio-alignment-sorting
description: Sort, merge, collate, and verify alignment-file order with samtools and pysam. Use when preparing BAM files for indexing, variant calling, or paired-end analysis.
tool_type: cli
primary_tool: samtools
license: MIT
author: GPTomics
---

## Version Compatibility

Reference examples tested with: pysam 0.22+, samtools 1.19+ (checked on samtools 1.24, pysam 0.24.1, Picard 3.5.0, bwa 0.7.19)

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Alignment Sorting

Sort alignment files by coordinate or read name using samtools and pysam.

## Reference Files

- Read [performance-and-errors.md](references/performance-and-errors.md) when choosing compression or diagnosing a failed sort. It contains the measured compression comparison and the Common Errors table.

**"Sort a BAM file"** -> Reorder reads by genomic coordinate (for indexing/variant calling) or by name (for paired-end processing).
- CLI: `samtools sort -o sorted.bam input.bam`
- Python: `pysam.sort('-o', 'sorted.bam', 'input.bam')`

## Sort Orders

| Order | Flag | Header written | Use Case |
|-------|------|----------------|----------|
| Coordinate | default | `SO:coordinate` | Indexing, visualization, variant calling |
| Name, natural | `-n` | `SO:queryname SS:queryname:natural` | samtools-only paired-end steps (fixmate, collate-style grouping) |
| Name, ASCII | `-N` | `SO:queryname SS:queryname:lexicographical` | Anything read by Picard, GATK or htsjdk |
| Tag | `-t TAG` | `SO:unsorted SS:unsorted:TAG:coordinate` | Group by a tag value (position is the secondary key) |
| Template-coordinate | `--template-coordinate` | `SO:unsorted SS:unsorted:template-coordinate GO:query` | fgbio GroupReadsByUmi (needs `MC` tags: run `samtools fixmate -m` first) |

Only coordinate-sorted files can be indexed. `-n`, `-N` and `-t` output fail `samtools index` (`hts_idx_push` errors).

### Name sort: `-n` vs `-N`

In samtools 1.24 `-n` is **natural** order (`read2` before `read10`), not lexicographic. `-N` is ASCII order (`read10` before `read2`), which is what Picard writes and expects for queryname. When those orders differ (as they do for most numeric-suffix read names), Picard 3.5.0 `MarkDuplicates` on `sort -n` output dies with `Alignments added out of order ... Sort order is queryname`, and `ValidateSamFile` reports `RECORD_OUT_OF_ORDER`; the same reads sorted with `-N` pass (planted-duplicate BAM: 50 pairs, 100 reads flagged). If the orders coincide, `-n` can be accepted, but that is not a portable assumption. Rule: use `-N` whenever a Picard/GATK/htsjdk tool will read the name-sorted file (or sort with Picard `SortSam`); `-n` is fine when only samtools reads it.

## samtools sort

### Sort by Coordinate (Default)
```bash
samtools sort -o sorted.bam input.bam
# Sort and index in one pass (writes sorted.bam.csi, not .bai)
samtools sort --write-index -o sorted.bam input.bam
```

### Sort by Read Name
```bash
samtools sort -n -o namesorted.bam input.bam    # natural; samtools-only consumers
samtools sort -N -o namesorted.bam input.bam    # ASCII; Picard/GATK/htsjdk consumers
```

### Multi-threaded Sorting
```bash
samtools sort -@ 8 -o sorted.bam input.bam
```
`-@ N` requests N additional worker threads. The number observed by a process monitor depends on the data, host, and samtools implementation, so do not treat it as a fixed `N+1` total. Output records are identical for any `-@` setting.

### Control Memory Usage
```bash
samtools sort -m 4G -@ 4 -o sorted.bam input.bam
```
`-m` is a per-thread budget, but observed peak memory also depends on data, spills, and implementation; measure it on your host instead of treating (`-@` + 1) x `-m` as exact. Values below `1M` are rejected. Data beyond the budget spills to temp files and is merged; the result is unchanged.

### Set Temporary Files
```bash
samtools sort -T /tmp/sort_tmp -o sorted.bam input.bam
```
`-T` is a **prefix**. For a non-existent path, temp files are named like `PREFIX.nnnn.bam`; if the supplied path is an existing directory, samtools creates its generated temporary files inside that directory. Put either form on fast local storage.

Slow sort: raise `-@` and `-m`, put `-T` on fast local disk, and consult [performance-and-errors.md](references/performance-and-errors.md) before choosing final-output compression.

### Specify Output Format
```bash
# Output as BAM (default)
samtools sort -O bam -o sorted.bam input.bam

# Output as CRAM
samtools sort -O cram --reference ref.fa -o sorted.cram input.bam
```

### Sort by Tag
```bash
# Sort by cell barcode (10x Genomics); output cannot be indexed
samtools sort -t CB -o sorted_by_barcode.bam input.bam
```

### Sort for fgbio (template-coordinate)
```bash
samtools fixmate -m namesorted.bam fixmate.bam     # adds the MC tag the sort needs
samtools sort --template-coordinate -o tc.bam fixmate.bam
```
Without `MC` tags the sort stops with `no MC tag. Please run samtools fixmate on file first.`

## samtools collate vs sort -n

| Tool | Algorithm | Output guarantee |
|------|-----------|------------------|
| `sort -n` / `sort -N` | Full sort by QNAME | Strict total order by name |
| `collate` | Hash-bucket grouping | Mates adjacent; between-mate order undefined (`SO:unsorted GO:query`) |

Use `collate` when a tool only needs mates adjacent (paired FASTQ extraction, re-aligning, fixmate, markdup pre-processing). It is not reliably faster: on a 192k-read slice on local WSL storage, `collate` took 0.74 s against 0.66 s for `sort -n` (samtools 1.24); on a network mount it was slower. Time both on the storage that will hold your data. `sort -n`/`-N` are only needed when a tool requires a real name order (see the table below).

```bash
# Paired FASTQ extraction
samtools collate -O -u in.bam tmp_prefix | \
    samtools fastq -1 R1.fq.gz -2 R2.fq.gz -0 /dev/null -s /dev/null -n -
```

### Sort Order Required by Downstream Tool

| Operation | Required sort |
|-----------|---------------|
| `samtools index` | coordinate (hard requirement) |
| `samtools fixmate -m` | name (`-n` or `-N`) or collate (needs mates adjacent); coordinate input stops with `Coordinate sorted, require grouped/sorted by queryname` |
| `samtools markdup` | coordinate (after fixmate) |
| Picard MarkDuplicates / ValidateSamFile | coordinate, or queryname in `-N` (ASCII) order; `-n` is rejected when natural and ASCII orders differ, so use `-N` |
| GATK MarkDuplicatesSpark | coordinate or queryname |
| `samtools mpileup` / `bcftools mpileup` | coordinate (name-sorted input stops with `The input is not sorted`) |
| GATK HaplotypeCaller / Mutect2 | coordinate and indexed |
| htseq-count | `-r pos` for coordinate-sorted, `-r name` for name-sorted; in a one-gene GTF test over the DNA test BAM, `-r name` on coordinate-sorted input counted 5599 vs 2820 reads with `-r pos`, so do not treat that ratio as universal (`-p` here is `--samout-format`, not "paired") |
| featureCounts † | coordinate or name; `-p` for paired-end |
| umi_tools dedup | coordinate (with index) |
| fgbio GroupReadsByUmi | any order accepted (template-coordinate recommended to avoid an internal re-sort; run `samtools fixmate -m`/`fgbio SetMateInformation` first) |
| fgbio CallMolecularConsensusReads | template-coordinate order, as emitted by `fgbio GroupReadsByUmi` after mate information is set |
| Sniffles, cuteSV, Manta, Delly † | coordinate and indexed |
| Salmon alignment-mode, RSEM (STAR `--quantMode TranscriptomeSAM`) † | mates adjacent, grouped by name; not verified whether a strict lexicographic order is needed, so `sort -N` is the safe choice |

† Tool not installed in the checking environment: the requirement is unverified here; confirm against the tool's own documentation.

## Check Sort Order

### From Header (a hint only)
```bash
samtools view -H input.bam | grep "^@HD"
# SO:coordinate = coordinate sorted
# SO:queryname = name sorted (SS: says natural or lexicographical)
# SO:unsorted = not sorted (tag and template-coordinate sorts also say unsorted)
```
A header can be wrong: a shuffled BAM labelled `SO:coordinate` passes the header check and then fails to index. Check the records.

### Verify Records
```bash
# Fails ("Unsorted positions on sequence ...") if reads are out of order within a contig
samtools index input.bam /tmp/check.bai && echo "indexable"
# After a long sort, catch truncated/empty output (rc 16 truncated, rc 4 zero-byte)
samtools quickcheck -v sorted.bam
```
`samtools index` accepts contigs in a different order from the `@SQ` header, so for a full check walk the records with pysam (below).

## pysam Python Alternative

### Sort with pysam
```python
import pysam

pysam.sort('-o', 'sorted.bam', 'input.bam')
pysam.sort('-n', '-o', 'namesorted.bam', 'input.bam')   # '-N' for Picard/GATK consumers
pysam.sort('-@', '4', '-m', '2G', '-o', 'sorted.bam', 'input.bam')
```

### Avoid In-Python Sorting

Do not load BAM records into a list and call `sorted()`. `pysam.sort()` calls samtools' external-merge sort which spills to disk; loading reads into memory blows up around ~30M reads (~10 GB human BAM). Always delegate to `pysam.sort()`:
```python
import pysam

pysam.sort('-@', '4', '-m', '2G', '-T', '/tmp/sortpfx',
           '-o', 'sorted.bam', 'input.bam')
```

### Check Sort Order in pysam
```python
import pysam

def is_coordinate_sorted(path, reference_filename=None):
    """True only if records are in (reference order, POS) order, with unmapped reads last."""
    with pysam.AlignmentFile(path, 'rb', check_sq=False,
                             reference_filename=reference_filename) as bam:
        if not bam.references:       # uBAM has no coordinate order to verify
            return False
        prev = (-1, -1)
        for r in bam.fetch(until_eof=True):
            key = (r.reference_id, r.reference_start) if r.reference_id >= 0 else (float('inf'), float('inf'))
            if key < prev:
                return False
            prev = key
    return True

def ensure_coordinate_sorted(input_bam, output_bam, reference_filename=None):
    if is_coordinate_sorted(input_bam, reference_filename):
        return input_bam
    pysam.sort('-o', output_bam, input_bam)
    return output_bam
```
This catches a mislabelled `SO:coordinate` BAM and reference blocks in the wrong order; reading `bam.header['HD']['SO']` alone does not. For CRAM, set `REF_PATH` or pass its reference as `reference_filename`; otherwise pysam cannot decode records and may report a misleading truncation error.

### Stream Sort from Aligner
For streaming from aligners, use shell pipes, and make the shell fail if any stage fails (`pipefail`); with `shell=True` alone a crashed aligner still returns 0 and leaves a truncated BAM:
```python
import subprocess

subprocess.run(
    ['bash', '-o', 'pipefail', '-c', 'bwa mem ref.fa reads.fq | samtools sort -o aligned.bam'],
    check=True
)
```

## samtools merge

Combine multiple BAM files into one. `samtools merge` does NOT validate sort-order consistency across inputs; mismatched inputs silently produce a malformed output (exit 0, header `SO:coordinate`, records out of order).

### Verify Sort Order Consistency First
```bash
for f in *.bam; do samtools view -H "$f" | head -1; done | sort -u
# Should print exactly ONE line, e.g. "@HD VN:1.6 SO:coordinate"
# (name-sorted: -n and -N inputs print different SS: lines; do not mix them)
```

### Merge (dedup identical @RG and @PG)
When merging BAMs from different lanes / machines / aligners, RG IDs may collide. `-c` deduplicates `@RG` records **by ID**, so two different read groups that share an ID (different SM/PU) collapse into one and their reads are mislabelled. Make colliding IDs unique upstream (`samtools addreplacerg`) before merge, otherwise GATK BQSR (which keys models by RGID/PU) silently produces wrong recalibration.
```bash
# -c deduplicates @RG records; -p deduplicates @PG records (samtools-merge(1))
samtools merge -c -p -@ 8 merged.bam sample1.bam sample2.bam sample3.bam
```

### Merge Name-sorted Inputs
```bash
# Pass the same flag the inputs were sorted with; a plain merge does not give name order
samtools merge -n -o merged_n.bam a_n.bam b_n.bam    # natural (sort -n inputs)
samtools merge -N -o merged_N.bam a_N.bam b_N.bam    # ASCII (sort -N inputs)
```

### Merge with Threads / from File List
```bash
samtools merge -@ 4 merged.bam sample1.bam sample2.bam sample3.bam
samtools merge -b files.txt merged.bam   # one BAM path per line
```

### Force Overwrite
```bash
samtools merge -f merged.bam sample1.bam sample2.bam
```

### Merge Specific Region
```bash
# -R needs every input indexed ("Could not retrieve index file" otherwise)
samtools merge -R chr1:1000000-2000000 merged_region.bam sample1.bam sample2.bam
```

### pysam Merge
```python
import pysam

pysam.merge('-c', '-p', '-f', 'merged.bam', 'sample1.bam', 'sample2.bam', 'sample3.bam')
```

## Common Workflows

**Goal:** Combine sorting with other alignment processing steps into efficient pipelines.

**Approach:** Pipe aligner output directly into `samtools sort` to avoid writing unsorted intermediates, then index for downstream access.

### Align and Sort
```bash
# One-time prerequisite: bwa index ref.fa
set -o pipefail    # without it a failed bwa in the pipe still exits 0
bwa mem -t 8 ref.fa R1.fq R2.fq | samtools sort -@ 4 -o aligned.bam
samtools index aligned.bam
```
`examples/sort_pipeline.sh` wraps this with input/index checks and a read-count check.

### Duplicate Marking
```bash
# collate groups mates (sort -n / sort -N work in its place); fixmate adds ms/MC; -u skips compression between steps
samtools collate -O -u in.bam tmp_prefix | \
    samtools fixmate -m -u - - | \
    samtools sort -u - | \
    samtools markdup - out.bam
samtools index out.bam
```

## Related Skills

- sam-bam-basics - View and convert alignment files
- alignment-indexing - Index after coordinate sorting
- duplicate-handling - Requires name-sorted input for fixmate
- alignment-filtering - Filter before or after sorting
