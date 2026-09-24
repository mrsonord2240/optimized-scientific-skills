## pysam Python Alternative

Two inputs break a plain `AlignmentFile(path, 'rb')`: an unaligned BAM has no `@SQ` lines (`ValueError: file has no sequences defined`), so the read-only snippets pass `check_sq=False`; and a CRAM needs `reference_filename='ref.fa'`, without which iterating fails with `OSError: truncated file` (region and pileup calls also need a `.crai`; `samtools flagstat` needs no reference, `samtools stats` needs `--reference`). `examples/qc_report.py` takes the reference as an optional second argument and exits 1 with a message when a CRAM cannot be decoded.

### Count Reads

**Goal:** Reproduce `samtools flagstat` (QC-passed column) counts and rates.

**Approach:** Count secondary/supplementary records and QC-failed reads separately; all rates use QC-passed primary reads only, so they equal flagstat's `primary`, `primary mapped`, `properly paired` and `primary duplicates` lines. `examples/qc_report.py` is the full version.

```python
import pysam

with pysam.AlignmentFile('input.bam', 'rb', check_sq=False) as bam:
    primary = mapped = paired = proper = 0
    for read in bam:
        if read.is_secondary or read.is_supplementary or read.is_qcfail:
            continue
        primary += 1
        mapped += not read.is_unmapped
        paired += read.is_paired
        proper += read.is_proper_pair

if not primary:
    raise SystemExit('no QC-passed primary reads')
print(f'Primary QC-passed: {primary}')
print(f'Mapped: {mapped} ({mapped/primary*100:.2f}% of primary)')
if paired:
    print(f'Properly paired: {proper} ({proper/paired*100:.2f}% of primary paired)')
else:
    print('Single-end data: no paired reads')
```

### Per-Chromosome Counts
```python
import pysam

with pysam.AlignmentFile('input.bam', 'rb') as bam:
    # needs an index (ValueError otherwise); counts include secondary/supplementary like idxstats
    for stat in bam.get_index_statistics():
        print(f'{stat.contig}: {stat.mapped} mapped, {stat.unmapped} unmapped')
```

### Depth in a Region (and at a Position)

**Goal:** Per-base depth, mean depth and breadth for a region that agree with `samtools depth -a`.

**Approach:** Fill a zero array over the region and record each pileup column; the mean and breadth divide by the region length, never by the number of covered columns.

**Reference (pysam 0.24.1, checked equal to `samtools depth -a` to 1e-9):**
```python
import pysam

def region_depth_stats(bam_path, chrom, start, end, thresholds=(10, 20), reference=None):
    """Depth over the 0-based half-open region [start, end); every base is in the denominator.
    A single position is the 1-bp region (pos - 1, pos). `reference` is the FASTA for a CRAM."""
    length = end - start
    if length <= 0:
        raise ValueError(f'empty region [{start}, {end})')
    depths = [0] * length
    with pysam.AlignmentFile(bam_path, 'rb', reference_filename=reference) as bam:
        # truncate=True: only columns inside the region (otherwise whole read footprints are returned)
        # max_depth: pysam's default cap is 8000; ignore_orphans=False: default drops paired reads
        #   that lack the proper-pair flag; min_base_quality=0: default is 13 (samtools depth: 0)
        # pileup.n counts deletions and ref-skips, so count real bases explicitly
        for col in bam.pileup(chrom, start, end, truncate=True, max_depth=1_000_000,
                              ignore_orphans=False, ignore_overlaps=False, min_base_quality=0):
            depths[col.reference_pos - start] = sum(1 for r in col.pileups if not r.is_del and not r.is_refskip)
    stats = {'length': length, 'covered': sum(d > 0 for d in depths),
             'mean_depth': sum(depths) / length, 'max_depth': max(depths)}
    stats['pct_covered'] = stats['covered'] / length * 100
    for t in thresholds:
        stats[f'pct_ge_{t}x'] = sum(d >= t for d in depths) / length * 100
    return stats

stats = region_depth_stats('input.bam', 'chr1', 1000000, 2000000)
print(f'Coverage: {stats["pct_covered"]:.1f}%  Mean depth: {stats["mean_depth"]:.1f}x  >=20x: {stats["pct_ge_20x"]:.1f}%')
```
Like `samtools depth`, this skips unmapped, secondary, QC-failed and duplicate reads and counts overlapping mates twice (the mean is 2x the `mosdepth` default on overlapping pairs). Python is slow beyond a few Mb: use `samtools coverage` or `mosdepth` for large regions.

### Insert Size Distribution

**Goal:** Compute the insert size distribution to assess library preparation quality.

**Approach:** Iterate QC-passed primary properly paired read1 records, accumulate template lengths into a Counter, then compute summary statistics.

**Reference (pysam 0.24.1):**
```python
import pysam
from collections import Counter

insert_sizes = Counter()

with pysam.AlignmentFile('input.bam', 'rb', check_sq=False) as bam:
    for read in bam:
        if read.is_secondary or read.is_supplementary or read.is_qcfail:
            continue
        if read.is_proper_pair and read.is_read1 and read.template_length > 0:
            insert_sizes[read.template_length] += 1

if not insert_sizes:
    raise SystemExit('no properly paired read1 records (single-end data, or proper-pair flag unset: see Insert Size Caveats)')
sizes = list(insert_sizes.keys())
mean_insert = sum(s * c for s, c in insert_sizes.items()) / sum(insert_sizes.values())
print(f'Mean insert size: {mean_insert:.0f}')
print(f'Min: {min(sizes)}, Max: {max(sizes)}')
```
