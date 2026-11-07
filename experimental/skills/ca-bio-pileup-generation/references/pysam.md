## pysam Python Alternative

`bam.pileup(contig, start, stop)` is 0-based, half-open; pass `truncate=True` or it also yields the columns of every overlapping read outside the region. `pileup_column.pos` is 0-based.

**Depth is `len(pileup_column.pileups)`, not `pileup_column.n`.** `n` counts reads before the base-quality filter and overlap removal, so it exceeds the mpileup depth (on the human test BAM it differs at 1087 of 1157 positions). `len(pileup_column.pileups)` (= `get_num_aligned()`) equals the `samtools mpileup` depth column when the parameters below match.

**BAQ is switched on by `fastafile`, and only under pysam's default stepper.** The default stepper is `'samtools'` (not `'all'`); leave `stepper` unset or pass `'samtools'`. Without `fastafile`, pysam defaults equal `samtools mpileup -B`. With `fastafile=pysam.FastaFile(ref)` BAQ is applied and the pileup equals `samtools mpileup -f ref.fa`; add `compute_baq=False` to get `-B` back. An explicit `stepper='all'` applies no BAQ even with `fastafile` (it equals `-B`) and also drops overlap and orphan handling, so it matches neither mpileup default. Checked position by position (depth and base counts) on 6 BAMs (human DNA, spliced RNA-seq, 1000G, ARTIC nanopore, 2 synthetic): 0 differing positions for each mapping in the table.

| samtools mpileup | `bam.pileup()` argument | Note |
|------------------|-------------------------|------|
| `-f ref.fa` (BAQ on) | `fastafile=FastaFile(ref)` | default stepper (`'samtools'`) only; `stepper='all'` gives no BAQ |
| `-B` | no `fastafile`, or `compute_baq=False` | |
| `-E` | `fastafile=FastaFile(ref), redo_baq=True` | recomputes over an existing BQ tag; without it the tag is reused |
| `-Q 13` (default) | `min_base_quality=13` (default) | pysam default matches |
| `-q N` | `min_mapping_quality=N` | |
| `-d N` (default 8000) | `max_depth=N` (default 8000) | `0` is not unlimited (see Maximum Depth) |
| `-x` | `ignore_overlaps=False` | default `True` matches |
| `-A` | `ignore_orphans=False` | default `True` matches |
| `--ff` | `flag_filter=INT` (default 1796 = the same four flags) | `--ff 0` = `flag_filter=0`; `stepper='nofilter'` also drops the orphan filter, so it is not `--ff 0` |
| `--rf` | `flag_require=INT` | any of the bits set, like `--rf` |
| `-C 50` | `adjust_capq_threshold=50` | |

Each read in `pileup_column.pileups` needs three checks in this order: **`is_refskip` first**, because pysam sets `is_del=True` on spliced-read `N` skips as well; then `is_del`; else the base at `query_position`. Testing `is_del` first reports intron positions as deletions. `pileup_read.indel` is the length of the insertion (>0) or deletion (<0) that follows the base.

### Basic Pileup
```python
import pysam

with pysam.AlignmentFile('input.bam', 'rb') as bam:
    for pileup_column in bam.pileup('chr1', 1000000, 1001000, truncate=True):
        print(f'{pileup_column.reference_name}:{pileup_column.pos + 1} depth={len(pileup_column.pileups)}')
```

### Access Reads at Position
```python
import pysam

with pysam.AlignmentFile('input.bam', 'rb') as bam:
    for pileup_column in bam.pileup('chr1', 1000000, 1000001, truncate=True):
        print(f'Position: {pileup_column.pos}')
        print(f'Depth: {len(pileup_column.pileups)}')

        for pileup_read in pileup_column.pileups:
            if pileup_read.is_refskip:
                print('  Reference skip')
            elif pileup_read.is_del:
                print('  Deletion')
            else:
                aln = pileup_read.alignment
                qpos = pileup_read.query_position
                strand = '-' if aln.is_reverse else '+'
                print(f'  {aln.query_name} {strand} {aln.query_sequence[qpos]} (Q{aln.query_qualities[qpos]})')
```

### Shipped pysam helpers

`scripts/pileup_helpers.py` contains the reusable `allele_counts`,
`allele_frequency`, `find_variants`, and `pileup_text` functions. For
`allele_counts`, use a **0-based pos: pass 1-based position minus 1**;
`find_variants` returns 1-based VCF-style positions. Read-base
`N`, reference skips, and reference-`N` sites are excluded from the
SNV-oriented counts, so frequency denominators agree across the helpers.

```python
from scripts.pileup_helpers import allele_counts, allele_frequency, find_variants

counts = allele_counts('input.bam', 'chr1', 1000000 - 1, min_base_quality=20)
frequency = allele_frequency('input.bam', 'chr1', 1000000 - 1, min_base_quality=20)
calls = find_variants('input.bam', 'reference.fa', 'chr1', 999000, 1000000)
```

Run the command-line counter from the Skill directory with
`python scripts/allele_counts.py input.bam chr1:1000000` (1-based, MAPQ/baseQ >= 20).


### Pileup with Quality Filtering
```python
import pysam

with pysam.AlignmentFile('input.bam', 'rb') as bam:
    for pileup_column in bam.pileup('chr1', 1000000, 1001000,
                                     truncate=True,
                                     min_mapping_quality=20,
                                     min_base_quality=20):
        print(f'{pileup_column.pos + 1}: {len(pileup_column.pileups)}')
```

### Generate Pileup Text

```python
from scripts.pileup_helpers import pileup_text

for row in pileup_text('input.bam', 'reference.fa', 'chr1', 1000000, 1000100):
    print(row)
```

`pileup_text` matches ordinary aligner CIGARs, including adjacent insertions
and deletions. It explicitly rejects padded (`P`) CIGARs; use `samtools
mpileup` when their text rendering is required. A deletion after a read's last
base is rendered with samtools' Q0 deletion quality rather than crashing.
