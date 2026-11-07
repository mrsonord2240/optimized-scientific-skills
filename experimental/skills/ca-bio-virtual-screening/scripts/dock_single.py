#!/usr/bin/env python3
"""Dock one ligand with AutoDock Vina and drop non-physical (positive-energy) modes.

Uses the Vina Python API (Vina 1.2+) when `from vina import Vina` works; otherwise falls back to
the Vina CLI through subprocess (Windows has no `vina` wheel: pip fails with "Boost library
location was not found"). Both paths take the same seed.

Inputs:  receptor PDBQT, ligand PDBQT, box centre and size (A), exhaustiveness, n poses, seed.
Output:  poses written to --out (PDBQT); prints one affinity per kept mode (kcal/mol).
Usage:   python dock_single.py receptor.pdbqt ligand.pdbqt --center X Y Z --size X Y Z \
             [--exhaustiveness 8] [--n-poses 10] [--seed 42] [--out poses.pdbqt] [--vina-exe vina]
"""
import argparse
import re
import subprocess


def dock_single(receptor_pdbqt, ligand_pdbqt, center, box_size, exhaustiveness=8,
                n_poses=10, seed=42, out_pdbqt='poses.pdbqt', vina_exe='vina'):
    try:
        from vina import Vina
    except ImportError:
        Vina = None
    if Vina is not None:
        v = Vina(sf_name='vina', seed=seed)
        v.set_receptor(receptor_pdbqt)
        v.set_ligand_from_file(ligand_pdbqt)
        v.compute_vina_maps(center=center, box_size=box_size)
        v.dock(exhaustiveness=exhaustiveness, n_poses=n_poses)
        v.write_poses(out_pdbqt, n_poses=n_poses, overwrite=True)
        energies = [list(e) for e in v.energies(n_poses=n_poses)]
    else:
        cmd = [vina_exe, '--receptor', receptor_pdbqt, '--ligand', ligand_pdbqt,
               '--center_x', str(center[0]), '--center_y', str(center[1]), '--center_z', str(center[2]),
               '--size_x', str(box_size[0]), '--size_y', str(box_size[1]), '--size_z', str(box_size[2]),
               '--exhaustiveness', str(exhaustiveness), '--num_modes', str(n_poses),
               '--seed', str(seed), '--out', out_pdbqt]
        run = subprocess.run(cmd, capture_output=True, text=True)
        if run.returncode != 0:
            raise RuntimeError(f'vina failed ({run.returncode}): {run.stderr or run.stdout}')
        energies = [[float(x) for x in m.groups()] for m in
                    re.finditer(r'REMARK VINA RESULT:\s+(\S+)\s+(\S+)\s+(\S+)', open(out_pdbqt).read())]
    # Vina occasionally emits a physically nonsensical positive-energy mode (e.g. +68 kcal/mol):
    # keep only affinity < 0.
    return [e for e in energies if e[0] < 0]


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('receptor_pdbqt')
    ap.add_argument('ligand_pdbqt')
    ap.add_argument('--center', type=float, nargs=3, required=True)
    ap.add_argument('--size', type=float, nargs=3, required=True)
    ap.add_argument('--exhaustiveness', type=int, default=8)
    ap.add_argument('--n-poses', type=int, default=10)
    ap.add_argument('--seed', type=int, default=42)
    ap.add_argument('--out', default='poses.pdbqt')
    ap.add_argument('--vina-exe', default='vina')
    a = ap.parse_args()
    for e in dock_single(a.receptor_pdbqt, a.ligand_pdbqt, a.center, a.size, a.exhaustiveness,
                         a.n_poses, a.seed, a.out, a.vina_exe):
        print(e[0])
