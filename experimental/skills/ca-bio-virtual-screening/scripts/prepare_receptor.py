#!/usr/bin/env python3
"""Prepare a receptor PDBQT: pdb2pqr protonation (PROPKA, AMBER) -> meeko.

Inputs:  a repaired PDB (waters/cofactors/metals already decided), target pH, output PDBQT path.
Output:  the receptor PDBQT (plus <base>_pH<pH>.pdb / .pqr beside the input).
Usage:   python prepare_receptor.py repaired.pdb receptor.pdbqt [--ph 7.4]
Needs:   pdb2pqr and mk_prepare_receptor (meeko's console script, no .py suffix) on PATH.
"""
import argparse
import subprocess
from pathlib import Path


def prepare_receptor(repaired_pdb, pdbqt_out, pH=7.4):
    # Decide which waters/cofactors/metals to retain before this function.
    base = str(Path(repaired_pdb).with_suffix(''))
    protonated_pdb = f'{base}_pH{pH}.pdb'
    pqr_file = f'{base}_pH{pH}.pqr'
    # --pdb-output writes a protonated PDB alongside the PQR; hand that PDB to
    # mk_prepare_receptor --read_pdb rather than --read_pqr (insertion-code residues
    # such as 184A crash meeko's PQR reader).
    subprocess.run(['pdb2pqr', '--ff=AMBER', f'--with-ph={pH}',
                    '--pdb-output', protonated_pdb,
                    repaired_pdb, pqr_file], check=True)
    output_basename = str(Path(pdbqt_out).with_suffix(''))
    subprocess.run(['mk_prepare_receptor', '--read_pdb', protonated_pdb,
                    '-o', output_basename, '-p'], check=True)
    return pdbqt_out


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('repaired_pdb')
    ap.add_argument('pdbqt_out')
    ap.add_argument('--ph', type=float, default=7.4)
    a = ap.parse_args()
    print(prepare_receptor(a.repaired_pdb, a.pdbqt_out, a.ph))
