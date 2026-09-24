'''Convert alignment between different formats'''
# Reference: biopython 1.83+ | Checked on biopython 1.88

import re
import sys
import tempfile
from pathlib import Path

from Bio import AlignIO

MOLECULE_TYPES = {'dna': 'DNA', 'rna': 'RNA', 'protein': 'protein'}
NUCLEOTIDE_CODES = frozenset('ACGTUNRYKMSWBDHVX')


def infer_molecule_type(alignment):
    '''Infer DNA, RNA, or protein from IUPAC residues; reject mixed T/U input.'''
    residues = ''.join(str(record.seq) for record in alignment).upper()
    residues = re.sub(r'[-.?*]', '', residues)
    if not residues:
        sys.exit('No residues in the alignment')
    if 'T' in residues and 'U' in residues:
        sys.exit('Mixed T and U residues: split or correct the alignment before choosing DNA or RNA')
    if not set(residues) <= NUCLEOTIDE_CODES:
        return 'protein'
    return 'RNA' if 'U' in residues else 'DNA'


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

    # NEXUS writes molecule_type into datatype= unchecked. Validate the complete
    # nucleotide alphabet and any override before producing any output files.
    inferred = infer_molecule_type(alignment)
    molecule_type = inferred
    if len(sys.argv) > 2:
        molecule_type = MOLECULE_TYPES.get(sys.argv[2].lower())
        if molecule_type is None:
            sys.exit(f"Unknown molecule type '{sys.argv[2]}': use DNA, RNA or protein")
        if molecule_type != inferred:
            sys.exit(f'Molecule type {molecule_type} contradicts the residues (look like {inferred})')
    print(f'Molecule type: {molecule_type}')

    unsafe_ids = [record.id for record in alignment if not re.fullmatch(r'[A-Za-z0-9_]+', record.id)]
    if unsafe_ids:
        print('Warning: NEXUS ids are not MrBayes-safe; sanitize characters outside [A-Za-z0-9_] before MrBayes.')

    for record in alignment:
        record.annotations['molecule_type'] = molecule_type

    # Write and validate a complete set in a temporary directory. This prevents a
    # rejected type or failed NEXUS write from leaving a mixture of fresh and empty files.
    with tempfile.TemporaryDirectory(dir='.') as temporary_dir:
        temporary_dir = Path(temporary_dir)
        for output_file, output_format in conversions:
            AlignIO.write(alignment, temporary_dir / output_file, output_format)

        expected = 'protein' if molecule_type == 'protein' else molecule_type.lower()
        nexus_text = (temporary_dir / 'output.nex').read_text()
        assert f'datatype={expected}' in nexus_text, f'output.nex is not datatype={expected}'
        reread = AlignIO.read(temporary_dir / 'output.nex', 'nexus')
        assert [str(record.seq) for record in reread] == [str(record.seq) for record in alignment]

        for output_file, output_format in conversions:
            (temporary_dir / output_file).replace(output_file)
            print(f'Wrote: {output_file} ({output_format})')

    print(f'Checked: output.nex says datatype={expected} and re-reads with identical sequences')
