#!/usr/bin/env python3
'''RDKit-only linker-reach screen -- a cheap, local pre-filter, NOT a ternary
complex prediction.

Implements SKILL.md's "Linker Geometry Assessment" `attachment_distance()` idea as a
runnable screen: for each candidate linker, sample 3D conformers and measure the
achievable distance between its two attachment points. Compare that range against a
required span (measured separately from your binary target-ligand / E3-ligand
co-crystal exit vectors) to flag which linkers can plausibly reach without assuming
a fully extended, unstrained geometry.

This has no knowledge of the target or E3 protein, binding-site sterics, or
interface energetics -- it only bounds linker reach. It does not replace PRosettaC,
DeepTernary, AlphaFold3, Boltz, or HADDOCK (see SKILL.md's "Ternary Complex
Prediction Tools" for how to run those externally); use it only to cut an
enumerated linker library down before spending an external submission on it.

Reference: RDKit 2024.09+ | Verify API if version differs
'''

import math
import sys
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem

sys.path.insert(0, str(Path(__file__).parent))
from protac_enumerate import LINKERS, _dummy_indices  # noqa: E402


def linker_reach(linker_smi, n_confs=30, seed=42):
    '''Sample 3D conformers of a standalone linker and measure the [*:1]-[*:2]
    attachment-point distance across them (min/mean/max, in Angstrom).
    '''
    mol = Chem.MolFromSmiles(linker_smi)
    if mol is None:
        raise ValueError('Linker SMILES is invalid')
    d1 = _dummy_indices(mol, atom_map=1)
    d2 = _dummy_indices(mol, atom_map=2)
    if len(d1) != 1 or len(d2) != 1:
        raise ValueError('Linker must contain exactly one [*:1] and one [*:2]')
    anchor1 = mol.GetAtomWithIdx(d1[0]).GetNeighbors()[0].GetIdx()
    anchor2 = mol.GetAtomWithIdx(d2[0]).GetNeighbors()[0].GetIdx()
    n_bonds = len(Chem.GetShortestPath(mol, d1[0], d2[0])) - 1

    removed = sorted([d1[0], d2[0]])
    new_anchor1 = anchor1 - sum(1 for r in removed if r < anchor1)
    new_anchor2 = anchor2 - sum(1 for r in removed if r < anchor2)

    rw = Chem.RWMol(mol)
    for idx in sorted(removed, reverse=True):
        rw.RemoveAtom(idx)
    capped = rw.GetMol()
    Chem.SanitizeMol(capped)
    capped = Chem.AddHs(capped)

    cids = AllChem.EmbedMultipleConfs(
        capped, numConfs=n_confs, randomSeed=seed,
        useRandomCoords=True, pruneRmsThresh=0.1,
    )
    if len(cids) == 0:
        raise RuntimeError(f'Conformer embedding failed for linker: {linker_smi}')
    AllChem.MMFFOptimizeMoleculeConfs(capped, maxIters=500)

    distances = [
        capped.GetConformer(cid).GetAtomPosition(new_anchor1).Distance(
            capped.GetConformer(cid).GetAtomPosition(new_anchor2)
        )
        for cid in cids
    ]
    return {
        'n_bonds': n_bonds,
        'min_A': min(distances),
        'mean_A': sum(distances) / len(distances),
        'max_A': max(distances),
        'n_confs': len(cids),
    }


def screen_linker_library(required_span_A, tolerance_A=2.0, linkers=None):
    '''Flag which named linkers (from protac_enumerate.LINKERS by default) have a
    sampled reach range overlapping `required_span_A +/- tolerance_A`.

    `required_span_A` should come from measuring `attachment_distance()` on your own
    binary co-crystal structures, not from this script.
    '''
    linkers = LINKERS if linkers is None else linkers
    results = {}
    for name, smi in linkers.items():
        reach = linker_reach(smi)
        feasible = (reach['min_A'] - tolerance_A) <= required_span_A <= (reach['max_A'] + tolerance_A)
        reach['feasible'] = feasible
        results[name] = reach
    return results


def _theoretical_max_reach_alkyl(n_bonds, bond_len=1.54, tetra_half_angle_deg=54.75):
    '''All-trans zigzag projection -- a classic extended-chain-length estimate, used
    below only as a sanity upper bound on the sampled reach for pure-alkyl linkers.
    '''
    return n_bonds * bond_len * math.sin(math.radians(tetra_half_angle_deg))


if __name__ == '__main__':
    # ILLUSTRATIVE required span -- replace with the real exit-vector distance
    # measured from your target-ligand / E3-ligand binary co-crystal structures.
    required_span_A = 9.0

    results = screen_linker_library(required_span_A)
    print(f'Required span (illustrative): {required_span_A} A\n')
    print(f'{"linker":<20}{"n_bonds":>8}{"min_A":>8}{"mean_A":>8}{"max_A":>8}{"feasible":>10}')
    for name, r in results.items():
        print(f'{name:<20}{r["n_bonds"]:>8}{r["min_A"]:>8.2f}{r["mean_A"]:>8.2f}'
              f'{r["max_A"]:>8.2f}{str(r["feasible"]):>10}')

    # Verification: reach must increase monotonically with chain length for the
    # pure-alkyl series, and must stay at or below the all-trans theoretical bound
    # (random-coordinate conformer sampling should not exceed full extension).
    alkyl_order = ['short_alkyl', 'medium_alkyl', 'long_alkyl']
    alkyl_max = [results[name]['max_A'] for name in alkyl_order]
    assert alkyl_max == sorted(alkyl_max), 'expected reach to grow with alkyl chain length'
    for name in alkyl_order:
        theo = _theoretical_max_reach_alkyl(results[name]['n_bonds'])
        assert results[name]['max_A'] <= theo + 0.5, (
            f'{name}: sampled max {results[name]["max_A"]:.2f} A exceeds the '
            f'all-trans theoretical bound {theo:.2f} A'
        )
    print('\nMonotonic-reach and theoretical-bound checks: PASS')
