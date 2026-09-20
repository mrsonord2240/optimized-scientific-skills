'''Convert alignment between different formats'''
# Reference: biopython 1.83+ | Checked on biopython 1.88

import sys
from pathlib import Path

from Bio import AlignIO

if __name__ == '__main__':
    input_file = sys.argv[1] if len(sys.argv) > 1 else str(Path(__file__).parent / 'sample_alignment.aln')
    input_format = 'clustal'
    molecule_type = sys.argv[2] if len(sys.argv) > 2 else 'DNA'  # 'DNA', 'RNA' or 'protein'; NEXUS output requires it

    conversions = [
        ('output.fasta', 'fasta'),
        ('output.phy', 'phylip-relaxed'),
        ('output.nex', 'nexus'),
    ]

    alignment = AlignIO.read(input_file, input_format)
    print(f'Read alignment: {len(alignment)} sequences, {alignment.get_alignment_length()} columns')

    for record in alignment:
        record.annotations['molecule_type'] = molecule_type

    for output_file, output_format in conversions:
        AlignIO.write(alignment, output_file, output_format)
        print(f'Wrote: {output_file} ({output_format})')
