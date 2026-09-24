"""Compare ETKDGv3's macrocycle default with an explicit opt-out.

Usage:
  python examples/compare_macrocycle_embedding.py --smiles 'O=C1CCCCCCCCCCCNC1'

This is a molecule-specific diagnostic, not a claim that either setting is
universally superior. It writes a JSON summary to stdout for audit or workflow
records.
"""

import argparse
import json

from rdkit import Chem
from rdkit.Chem import AllChem


def embed_and_score(smiles: str, n_conf: int, seed: int, macrocycle_torsions: bool) -> dict:
    """Embed and MMFF94s-score one ETKDGv3 setting deterministically."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None or mol.GetNumAtoms() == 0:
        raise ValueError(f"Invalid SMILES: {smiles!r}")
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    params.useRandomCoords = True
    params.maxIterations = 5000
    params.useMacrocycleTorsions = macrocycle_torsions
    ids = list(AllChem.EmbedMultipleConfs(mol, numConfs=n_conf, params=params))

    properties = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant="MMFF94s")
    if properties is None:
        raise ValueError("MMFF94s parameters are unavailable for this comparison molecule")
    energies = []
    for conf_id in ids:
        force_field = AllChem.MMFFGetMoleculeForceField(mol, properties, confId=conf_id)
        if force_field is None:
            continue
        force_field.Minimize(maxIts=1000)
        energies.append(float(force_field.CalcEnergy()))
    return {
        "requested": n_conf,
        "embedded": len(ids),
        "mmff_scored": len(energies),
        "energy_range_kcal_mol": None if not energies else [min(energies), max(energies)],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--smiles", required=True, help="Macrocycle SMILES to compare")
    parser.add_argument("--n-conf", type=int, default=40, help="Requested conformers per setting")
    parser.add_argument("--seed", type=int, default=42, help="RDKit random seed")
    args = parser.parse_args()
    if args.n_conf < 1:
        raise ValueError("--n-conf must be at least 1")

    default = embed_and_score(args.smiles, args.n_conf, args.seed, True)
    opt_out = embed_and_score(args.smiles, args.n_conf, args.seed, False)
    print(json.dumps({"smiles": args.smiles, "seed": args.seed,
                      "etkdgv3_default": default,
                      "macrocycle_torsions_opt_out": opt_out}, indent=2))


if __name__ == "__main__":
    main()
