# MAF Block Coordinate Conventions

Read this reference when converting UCSC MAF coordinates to plus-strand genome coordinates.

UCSC MAF (read via `AlignIO.parse(file, 'maf')`) returns blocks with per-row `annotations`:

- `start` (0-based; converts directly to BED but is off-by-one vs GFF)
- `size` (length on src strand)
- `strand` (an **int**, `1` or `-1`, not the file's `+` / `-` characters; comparing to `'-'` never matches and silently returns unconverted minus-strand starts)
- `srcSize` (length of source chromosome)

For minus-strand rows, `start` is measured from the END of the source contig: the corresponding plus-strand start is `srcSize - start - size`. Without this conversion, lifting MAF to genome coordinates places minus-strand blocks at the wrong locus. Reference: UCSC MAF spec at genome.ucsc.edu/FAQ/FAQformat.html#format5.

```python
def maf_to_plus_strand_coords(row_anno):
    if row_anno['strand'] == -1:
        return row_anno['srcSize'] - row_anno['start'] - row_anno['size']
    return row_anno['start']

# Ground-truth check against your own reference (contigs: dict name -> str). Bio.Seq.reverse_complement
# of the plus-strand slice must equal the ungapped minus-strand row.
for block in AlignIO.parse('blocks.maf', 'maf'):
    for record in block:
        a = record.annotations
        start = maf_to_plus_strand_coords(a)
        fragment = Seq(contigs[record.id.split('.', 1)[1]][start:start + a['size']])
        if a['strand'] == -1:
            fragment = fragment.reverse_complement()
        assert str(fragment).upper() == str(record.seq).replace('-', '').upper()
```
