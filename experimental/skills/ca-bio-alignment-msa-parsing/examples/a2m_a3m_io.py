'''Read A2M / A3M alignments (HMMER, HHsuite, ColabFold) and extract match-only columns.

A2M: each sequence has identical match-column count; lowercase = insert state, '.' = gap
in an insert column, '-' = gap in a match column. HH-suite writes A2M padded to a rectangle, but
HMMER 3.4 `hmmalign --outformat a2m` writes UNPADDED rows (no '.' characters, unequal lengths), which
AlignIO.read(..., 'fasta') rejects. This script therefore reads with SeqIO.parse, which accepts both.
A3M: insert columns are not padded; lowercase characters appear inline per sequence.
'''
# Reference: biopython 1.83+ | Verify API if version differs

from Bio import SeqIO

from msa_utils import example_path

def match_only_columns(a2m_records):
    # a2m_records: any iterable of SeqRecords (padded or unpadded rows).
    # Keep upper-case residues and '-' (match columns); drops lowercase inserts and their '.' padding.
    return [
        ''.join(c for c in str(record.seq) if c.isupper() or c == '-')
        for record in a2m_records
    ]

if __name__ == '__main__':
    alignment = list(SeqIO.parse(example_path('example.a2m'), 'fasta'))
    widths = sorted({len(record.seq) for record in alignment})
    print(f'A2M file: {len(alignment)} sequences, row lengths {widths} ({"padded" if len(widths) == 1 else "unpadded"})')

    match_columns = match_only_columns(alignment)
    counts = sorted({len(row) for row in match_columns})
    print(f'After dropping insert states: {counts} match columns per row')
    if len(counts) != 1:
        print('WARNING: rows differ in match-column count; this is not a valid A2M')
    for record, match_only in list(zip(alignment, match_columns))[:5]:
        inserts = sum(1 for c in str(record.seq) if c.islower())
        print(f'{record.id}: {len(match_only)} match, {inserts} inserts')
