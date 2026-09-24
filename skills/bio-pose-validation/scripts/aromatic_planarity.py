# Purpose: supplementary aromatic-ring planarity diagnostic (max out-of-plane deviation, Angstrom).
#          NOT a gate: PoseBusters' own `aromatic_ring_flatness` column is authoritative (see SKILL.md).
# Inputs:  an SDF file (one or more poses).
# Usage:   python scripts/aromatic_planarity.py poses.sdf
import sys

import numpy as np
from rdkit import Chem


def aromatic_planarity(mol):
    deviations = []
    for ring in mol.GetRingInfo().AtomRings():
        ring_atoms = [mol.GetAtomWithIdx(i) for i in ring]
        if not all(a.GetIsAromatic() for a in ring_atoms):
            continue
        coords = np.array([mol.GetConformer().GetAtomPosition(i)
                          for i in ring])
        centroid = coords.mean(axis=0)
        centered = coords - centroid
        _, s, vh = np.linalg.svd(centered)
        normal = vh[-1]
        deviation = np.abs(centered @ normal).max()
        deviations.append(deviation)
    return max(deviations) if deviations else 0


if __name__ == '__main__':
    sdf = sys.argv[1]
    for i, mol in enumerate(Chem.SDMolSupplier(sdf, removeHs=False)):
        if mol is None:
            print(f'pose {i}: parse failed')
            continue
        print(f'pose {i}: max aromatic ring deviation {aromatic_planarity(mol):.4f} A')
