'''Convert a Spectronaut library TSV to DIA-NN column names.

Conversion is renaming columns and reconciling units, not a copy: iRT and RelativeIntensity are
renamed, not recomputed, and the script refuses input whose iRT column is outside -50..200
(RT not in iRT units; check the column before converting).
Input:  Spectronaut library TSV with the columns in SPECTRONAUT_TO_DIANN.
Output: DIA-NN-named TSV.
Usage:  python spectronaut_to_diann.py --in spectronaut.tsv --out diann.tsv
'''
import argparse

import pandas as pd

# Spectronaut -> DIA-NN column mapping; iRT and RelativeIntensity are renamed, not recomputed.
SPECTRONAUT_TO_DIANN = {'ModifiedPeptide': 'ModifiedPeptide', 'iRT': 'iRT',
                        'RelativeIntensity': 'LibraryIntensity', 'FragmentMz': 'ProductMz',
                        'FragmentNumber': 'FragmentSeriesNumber', 'PrecursorMz': 'PrecursorMz',
                        'PrecursorCharge': 'PrecursorCharge', 'FragmentCharge': 'FragmentCharge',
                        'FragmentType': 'FragmentType', 'Genes': 'Genes'}


def spectronaut_to_diann(lib):
    out = lib.rename(columns=SPECTRONAUT_TO_DIANN)
    assert out['iRT'].between(-50, 200).all(), 'RT not in iRT units; check column before converting'
    return out


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--in', dest='inp', required=True)
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    lib = pd.read_csv(a.inp, sep='\t')
    out = spectronaut_to_diann(lib)
    out.to_csv(a.out, sep='\t', index=False)
    print(f'{len(out)} rows -> {a.out}')
