'''Convert alignment between different formats'''
# Reference: biopython 1.83+ | Checked on biopython 1.88

import re
import sys
from pathlib import Path

from Bio import AlignIO

MOLECULE_TYPES = {'dna': 'DNA', 'rna': 'RNA', 'protein': 'protein'}


def infer_molecule_type(alignment):
    '''Guess DNA / RNA / protein from the residues: >= 90% ACGTUN means nucleotide'''
    residues = ''.join(str(record.seq) for record in alignment).upper()
    residues = re.sub(r'[-.?*]', '', residues)
    if not residues:
        sys.exit('No residues in the alignment')
    if sum(residues.count(c) for c in 'ACGTUN') / len(residues) < 0.9:
        return 'protein'
    return 'RNA' if 'U' in residues and 'T' not in residues else 'DNA'


if __name__ == '__main__':
    input_file = sys.argv[1] if len(sys.argv) > 1 else str(Path(__file__).parent / 'sample_alignment.aln')
    input_format = 'clustal'

    conversions = [
        ('output.fasta', 'fasta'),
        ('output.phy', 'phylip-relaxed'),
        ('output.nex', 'nexus'),
    ]

    alignment = AlignIO.read(input_file, input_format)
    print(f'Read alignment: {len(alignment)} sequences, {alignment.get_alignment_length()} columns')

    # NEXUS needs a molecule type and writes it into `datatype=` unchecked, so infer it and
    # accept an optional override (argv[2]) only when it agrees with the residues
    inferred = infer_molecule_type(alignment)
    molecule_type = inferred
    if len(sys.argv) > 2:
        molecule_type = MOLECULE_TYPES.get(sys.argv[2].lower())
        if molecule_type is None:
            sys.exit(f"Unknown molecule type '{sys.argv[2]}': use DNA, RNA or protein")
        if (molecule_type == 'protein') != (inferred == 'protein'):
            sys.exit(f'Molecule type {molecule_type} contradicts the residues (look like {inferred})')
    print(f'Molecule type: {molecule_type}')

    for record in alignment:
        record.annotations['molecule_type'] = molecule_type

    for output_file, output_format in conversions:
        AlignIO.write(alignment, output_file, output_format)
        print(f'Wrote: {output_file} ({output_format})')

    # Check the NEXUS datatype line matches the alphabet
    expected = 'protein' if molecule_type == 'protein' else molecule_type.lower()
    nexus_text = Path('output.nex').read_text()
    assert f'datatype={expected}' in nexus_text, f'output.nex is not datatype={expected}'
    print(f'Checked: output.nex says datatype={expected}')
