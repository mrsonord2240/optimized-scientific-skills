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
| `wc` | Write CRAM (give `reference_filename=`; without it pysam warns and writes an embedded-reference CRAM) |

`bam.mapped` / `bam.unmapped` come from the BAM index and are unavailable for SAM, unindexed BAM and CRAM (0 or an error); count with a scan instead, as `examples/view_bam.py <file> [limit] [reference.fa]` does. On an unindexed CRAM, pysam prints `[E::cram_index_load]` lines on stderr; they are harmless when the run succeeds.

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
