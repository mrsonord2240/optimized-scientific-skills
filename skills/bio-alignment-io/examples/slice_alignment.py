'''Slice and subset alignments'''
# Reference: biopython 1.83+ | Checked on biopython 1.88

import sys
from pathlib import Path

from Bio import AlignIO

if __name__ == '__main__':
    input_file = sys.argv[1] if len(sys.argv) > 1 else str(Path(__file__).parent / 'sample_alignment.aln')
    col_start, col_end = 5, 15  # 0-based, end-exclusive

    alignment = AlignIO.read(input_file, 'clustal')
    n_cols = alignment.get_alignment_length()
    print(f'Original: {len(alignment)} sequences, {n_cols} columns')

    # Out-of-range column slices silently return 0 columns, so refuse them
    if col_end > n_cols:
        sys.exit(f'Columns {col_start}-{col_end} are outside the {n_cols}-column alignment')

    subset_seqs = alignment[0:5]
    print(f'First 5 sequences: {len(subset_seqs)} sequences')

    trimmed = alignment[:, col_start:col_end]
    print(f'Columns {col_start}-{col_end}: {trimmed.get_alignment_length()} columns')

    region = alignment[0:5, col_start:col_end]
    print(f'Combined slice: {len(region)} sequences, {region.get_alignment_length()} columns')

    AlignIO.write(region, 'trimmed_subset.fasta', 'fasta')
    print('Wrote trimmed subset to trimmed_subset.fasta')
