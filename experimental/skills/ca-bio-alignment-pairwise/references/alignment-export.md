# Substitution Counts and Export Formats (reference)

## Substitution Matrix from Alignment

**Goal:** Extract observed substitution frequencies from a completed alignment.

**Approach:** Access the `.substitutions` property to get a matrix of observed base/residue substitution counts.

```python
alignment = alignments[0]
substitutions = alignment.substitutions

# View as array (rows=target, cols=query)
print(substitutions)

# Access specific substitution counts
# substitutions['A', 'T'] gives count of A aligned to T
```

## Export Alignment to Different Formats

**Goal:** Convert an alignment to standard bioinformatics file formats for downstream tools.

**Approach:** Use Python's `format()` function with format specifiers (fasta, clustal, psl, sam) on the alignment object.

```python
alignment = alignments[0]

# Various output formats
print(format(alignment, 'fasta'))     # FASTA format
print(format(alignment, 'clustal'))   # Clustal format
print(format(alignment, 'psl'))       # PSL format (BLAT)
print(format(alignment, 'sam'))       # SAM format
```
