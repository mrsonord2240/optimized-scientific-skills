'''Build an OpenSWATH transition TSV (y-ion transitions with real fragment m/z) from a peptide list.

OpenSwathDecoyGenerator needs ALL of: a literal Annotation column, chemically real theoretical
fragment m/z, and a per-precursor grouping column (transition_group_id, unique per
PeptideSequence + PrecursorCharge). Without the grouping column TargetedFileConverter silently
merges peptides that share a charge. This script writes all three.
Input:  TSV with columns sequence, charge, protein, irt (header row required).
Output: OpenSWATH TSV; feed it to TargetedFileConverter (see references/format-conversion.md).
Usage:  python build_openswath_tsv.py --peptides peptides.tsv --out library.tsv [--n-frag 6]
Checked with pyteomics 5.0.1, OpenMS 3.5.0.
'''
import argparse

import pandas as pd
from pyteomics import mass


def build_openswath_tsv(peptides, path, n_frag=6):
    """peptides: [(sequence, charge, protein, iRT)] -> OpenSWATH TSV of y-ion transitions."""
    rows = []
    for seq, z, prot, irt in peptides:
        group = f'{seq}_{z}'  # unique per PeptideSequence + PrecursorCharge; required
        for i in range(1, n_frag + 1):
            rows.append({'PrecursorMz': mass.fast_mass(seq, charge=z),
                         'ProductMz': mass.fast_mass(seq[-i:], ion_type='y', charge=1),  # real m/z
                         'Tr_recalibrated': irt, 'transition_name': f'{group}_y{i}',
                         'transition_group_id': group, 'decoy': 0, 'LibraryIntensity': 1000.0 / i,
                         'PeptideSequence': seq, 'FullUniModPeptideName': seq,
                         'ProteinName': prot, 'PrecursorCharge': z, 'FragmentType': 'y',
                         'FragmentSeriesNumber': i, 'FragmentCharge': 1,
                         'Annotation': f'y{i}^1'})  # literal Annotation; required
    pd.DataFrame(rows).to_csv(path, sep='\t', index=False)


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--peptides', required=True, help='TSV: sequence, charge, protein, irt')
    ap.add_argument('--out', required=True)
    ap.add_argument('--n-frag', type=int, default=6)
    a = ap.parse_args()
    df = pd.read_csv(a.peptides, sep='\t')
    peptides = list(df[['sequence', 'charge', 'protein', 'irt']].itertuples(index=False, name=None))
    build_openswath_tsv(peptides, a.out, a.n_frag)
    print(f'{len(peptides)} precursors -> {a.out}')
