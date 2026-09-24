# Annotations and Position Mapping

Read when the request involves Stockholm annotations (SS_cons, RF, GS/GR) or converting between alignment columns and ungapped residue or PDB positions.

## Working with Annotations

Stockholm-derived alignments expose secondary-structure markup, GC/GR per-column annotations, and per-sequence metadata via `record.annotations`, `record.letter_annotations`, and `alignment.column_annotations`:

```python
alignment = AlignIO.read('pfam.sto', 'stockholm')
for record in alignment:
    if 'secondary_structure' in record.letter_annotations:
        print(record.id, record.letter_annotations['secondary_structure'])

ss_cons = alignment.column_annotations.get('secondary_structure')
```

GC SS_cons (consensus secondary structure), GC RF (reference coordinates), and GS metadata (organism, taxonomy) survive read/write through the `'stockholm'` format string but are silently discarded when writing to FASTA, PHYLIP, or NEXUS. Biopython 1.88 keeps only the tags it recognises even on a Stockholm round trip: on the Pfam PF00042 seed `#=GC seq_cons` and `#=GR pAS` were dropped (`#=GR AS` and `#=GS AC` kept, a `#=GS DE` line added), so keep the original file as the master copy if annotations matter for downstream analysis. `select_columns`, `remove_gappy_columns` and the filters above keep them; hand-built `SeqRecord(Seq(new_seq), id=...)` loops do not.

## Position Mapping

**Goal:** Convert between alignment column coordinates and ungapped sequence coordinates.

**Approach:** For one-off lookups, walk the sequence tracking gap characters. For repeated queries on the same sequence, vectorize with `numpy.cumsum` over a gap-mask -- O(L) preprocessing, O(1) lookups.

### Vectorized Coordinate Mapping (Recommended)

```python
import numpy as np

def coordinate_map(record):
    chars = np.frombuffer(str(record.seq).encode('ascii'), dtype=np.uint8)
    is_residue = ~np.isin(chars, [ord('-'), ord('.')])
    seq_to_aln = np.flatnonzero(is_residue)
    aln_to_seq = np.where(is_residue, np.cumsum(is_residue) - 1, -1)
    return seq_to_aln, aln_to_seq

seq_to_aln, aln_to_seq = coordinate_map(alignment[0])
column_of_residue_index_42 = seq_to_aln[42]   # 0-based index 42 is the 43rd residue
residue_index_at_column_100 = aln_to_seq[100]
```

`aln_to_seq[i] == -1` indicates a gap at alignment column `i`. This pattern handles 1 M-site genomic alignments in milliseconds compared to the loop-based version's seconds. It subsumes the single-lookup walk (`seq_pos += 1` for each non-gap char) in both directions.

### Mapping Alignment Columns to PDB Residues

A column-to-PDB mapping requires THREE coordinate systems: alignment column -> SEQRES residue (ungapped FASTA) -> ATOM residue (resolved structure). The SEQRES-to-ATOM map is non-trivial because PDB structures have unmodelled loops, N-terminal tags, engineered mutations, and seleno-substitutions. Conservation scores mapped via the bare alignment-to-SEQRES path will be off-by-many residues whenever the structure has missing density. Even the simplest case is offset: PDB 1MBN numbers residues from the first Val, so His93 is UniProt P02185 residue 94 (0-based index 93, checked). For SEQRES/ATOM extraction, `missing_residues`, and mapping by residue number (auth_seq_id) or SIFTS rather than string index, see `structural-biology/structure-navigation` (sections "Reading the Declared (SEQRES) Sequence and Locating Gaps" and "Reading mmCIF with an Explicit Numbering Scheme").
