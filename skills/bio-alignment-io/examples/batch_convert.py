'''Batch convert alignment files'''
# Reference: biopython 1.83+ | Checked on biopython 1.88

import sys
from pathlib import Path
from Bio import AlignIO

if __name__ == '__main__':
    input_dir = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).parent  # holds sample_alignment.aln
    output_dir = Path('converted/')

    input_format = 'clustal'
    output_format = 'fasta'

    input_files = sorted(input_dir.glob('*.aln'))
    if not input_files:
        sys.exit(f'No *.aln files found in {input_dir}')
    output_dir.mkdir(exist_ok=True)

    for input_file in input_files:
        alignment = AlignIO.read(input_file, input_format)
        output_file = output_dir / f'{input_file.stem}.fasta'
        AlignIO.write(alignment, output_file, output_format)
        print(f'{input_file.name} -> {output_file.name} ({len(alignment)} seqs)')

    print(f'\nConverted {len(input_files)} files to {output_format} format')
