## pysam Python Alternative

`examples/filter_bam.py` wraps the filters below as a command line (`python filter_bam.py in.bam out.bam -q 30 -d -p -P -r chr1:1000-2000`): region in samtools syntax (1-based, inclusive; contig, contig:start and commas accepted), warns when `-d` finds no duplicate flags, and indexes the output when it is coordinate-sorted.

### Filter with Function

**Goal:** Apply a multi-criteria quality filter to produce clean alignments for downstream analysis.

**Approach:** A predicate checking mapped status, primary alignment, duplicate flag, and MAPQ, streamed over the reads (`examples/filter_bam.py`). Equal to `samtools view -F 3332 -q 30` (record-for-record identical on two real BAMs).

`python examples/filter_bam.py input.bam filtered.bam -q 30 -d -p` (mapped, primary, non-duplicate, MAPQ >= 30; duplicates must already be marked).

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

**Goal:** Extract only reads overlapping target regions defined in a BED file, each read once, in coordinate order (same records as `samtools view -L` for tab- or space-delimited rows with start < end; `-L` reads a zero-width row, start == end, as "reads spanning that position strictly inside", which `fetch` cannot express, so drop or widen such rows).

**Approach:** Parse BED (skipping header lines), sort and merge the intervals, fetch each merged interval, and skip reads already written through the previous interval: such a read starts before the previous interval's end. A plain loop over BED rows writes a read once per overlapped row (802 duplicated records out of 3410 on the test BAM).

**Reference (pysam 0.22+):**
```bash
python scripts/filter_by_bed.py input.bam targets.bed targets.bam
```

### Subsample (Pair-Consistent)

Hash on QNAME so mates stay together (a fresh `random.random()` per read drops mates inconsistently and breaks paired-end tools). Mix the seed into a real hash: `crc32(qname) ^ seed` only flips the low bits and keeps the same reads for every seed, and `crc32(f'{seed}:{qname}')` still correlates between seeds (two seeds shared 65% of their reads). This picks different reads than `samtools view -s`.
```bash
python scripts/subsample_pysam.py input.bam subset.bam 0.1 42   # fraction, seed
```
