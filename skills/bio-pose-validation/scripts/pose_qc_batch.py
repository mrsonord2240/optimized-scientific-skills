# Purpose: PoseBusters `dock` QC over many docked-pose SDFs; keeps the first PB-valid pose per file.
# Inputs:  receptor PDB (with hydrogens) and one or more docked-pose SDFs, in the order the docker ranked them.
# Usage:   python scripts/pose_qc_batch.py receptor.pdb poses1.sdf [poses2.sdf ...]
import sys

import pandas as pd
from posebusters import PoseBusters


def pose_qc_pipeline(docked_sdfs, receptor_pdb):
    bust = PoseBusters(config='dock')
    all_results = []
    for sdf in docked_sdfs:
        r = bust.bust(mol_pred=sdf, mol_cond=receptor_pdb)
        check_cols = [
            col for col in r.select_dtypes(include='bool').columns
            if not col.lower().startswith('rmsd')
        ]
        r['pb_valid'] = r[check_cols].all(axis=1)
        r['source'] = sdf
        all_results.append(r)
    df = pd.concat(all_results)

    df['rank'] = df.groupby('source')['pb_valid'].cumsum()
    valid_top = df[df['pb_valid']].groupby('source').head(1)
    return valid_top


if __name__ == '__main__':
    receptor, sdfs = sys.argv[1], sys.argv[2:]
    top = pose_qc_pipeline(sdfs, receptor)
    print(f'{len(top)} / {len(set(sdfs))} files have a PB-valid pose')
    print(top[['source', 'rank', 'pb_valid']])
