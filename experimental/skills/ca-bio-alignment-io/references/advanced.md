# Advanced I/O



**Goal:** Build an alignment object from sequences defined in code rather than read from a file.

**Approach:** Construct SeqRecord objects with gap characters and wrap them in a MultipleSeqAlignment.

```python
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

records = [
    SeqRecord(Seq('ACTGACTGACTG'), id='seq1'),
    SeqRecord(Seq('ACTGACT-ACTG'), id='seq2'),
    SeqRecord(Seq('ACTG-CTGACTG'), id='seq3'),
]
alignment = MultipleSeqAlignment(records)
AlignIO.write(alignment, 'new_alignment.fasta', 'fasta')
```




**Goal:** Convert a directory of alignment files from one format to another in bulk.

**Approach:** Glob for input files and iterate, reading each alignment and writing to the target format.

```python
from pathlib import Path

input_dir = Path('alignments/')
output_dir = Path('converted/')
output_dir.mkdir(exist_ok=True)

for input_file in input_dir.glob('*.aln'):
    alignment = AlignIO.read(input_file, 'clustal')
    output_file = output_dir / f'{input_file.stem}.fasta'
    AlignIO.write(alignment, output_file, 'fasta')
```

## Alternative: Bio.Align Module I/O

**Goal:** Use the newer Bio.Align module for alignment I/O with access to features like counts and substitutions.

**Approach:** Use `Align.read()`, `Align.parse()`, and `Align.write()` which return `Alignment` objects instead of `MultipleSeqAlignment`. Stick to FASTA, Clustal, PHYLIP and MAF here; on Biopython 1.88 Stockholm through `Bio.Align` fails (see Common Errors) and its NEXUS writer needs `molecule_type`, so use `Bio.AlignIO` for those.

```python
from Bio import Align

# Read single alignment (returns Alignment object)
alignment = Align.read('alignment.aln', 'clustal')
print(alignment.shape, alignment.counts())   # (n_seqs, n_columns), AlignmentCounts

# Parse multiple alignments (e.g. MAF blocks)
for block in Align.parse('blocks.maf', 'maf'):
    print(f'Alignment with {len(block)} sequences')

# Write alignment
Align.write(alignment, 'output.fasta', 'fasta')
```

### When to Use Which

| Use Case | Module |
|----------|--------|
| Legacy code, MultipleSeqAlignment needed | `Bio.AlignIO` |
| Stockholm (annotations), NEXUS, `molecule_type` handling | `Bio.AlignIO` |
| Modern features (counts, substitutions) on FASTA / Clustal / PHYLIP / MAF | `Bio.Align` |
| Format conversion | `Bio.AlignIO` (works for every format above) |
| Working with pairwise alignments | `Bio.Align` |
