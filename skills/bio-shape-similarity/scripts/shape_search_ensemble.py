# Conformer-ensemble shape search (Open3DAlign) and scaffold-hop filter.
# Purpose: score each library molecule by the best O3A-aligned shape Tanimoto over its MMFF-converged
#          conformers against the query, report every dropped molecule, optionally keep only
#          shape-high / ECFP4-low scaffold-hop candidates.
# Inputs:  --query SMILES (connected), --library FILE.smi (one SMILES per line, optional name after a space)
# Usage:   python shape_search_ensemble.py --query 'CC(=O)Nc1ccc(C(=O)c2ccccc2)cc1' --library lib.smi [--n-conf 20]
#          python shape_search_ensemble.py --query ... --library lib.smi --scaffold-hop --shape-threshold 0.7 --ecfp-threshold 0.5
# Import:  from shape_search_ensemble import shape_search_ensemble, scaffold_hop_candidates
# Thresholds are repository starting defaults only; calibrate on a task-relevant active/decoy benchmark.
# Checked on RDKit 2026.03.6.
import argparse

from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, rdFingerprintGenerator, rdMolAlign, rdShapeHelpers


def shape_search_ensemble(query_mol, library_mols, n_conf=20, seed=42):
    hits = []
    dropped = []  # (smiles, reason) for every molecule that produced no usable score
    for target in library_mols:
        smi = Chem.MolToSmiles(target)
        if len(Chem.GetMolFrags(target)) > 1:
            dropped.append((smi, 'disconnected fragments (salt/multi-component); '
                                  'no single shape to compare'))
            continue

        target = Chem.AddHs(target)
        params = AllChem.ETKDGv3()
        params.randomSeed = seed
        ids = list(AllChem.EmbedMultipleConfs(target, numConfs=n_conf, params=params))
        if not ids:
            dropped.append((smi, 'embedding failed for all requested conformers'))
            continue
        if not AllChem.MMFFHasAllMoleculeParams(target):
            dropped.append((smi, 'MMFF parameters unavailable'))
            continue

        optimization = AllChem.MMFFOptimizeMoleculeConfs(target)
        # Keep only the conformer ids that converged (status == 0); a molecule with
        # 19 good conformers and 1 non-convergent one must not be discarded outright.
        converged_ids = [cid for cid, (status, _) in zip(ids, optimization) if status == 0]
        n_failed = len(ids) - len(converged_ids)
        if n_failed:
            print(f'WARNING: {smi}: {n_failed}/{len(ids)} conformers failed MMFF '
                  f'convergence; scoring the remaining {len(converged_ids)}')
        if not converged_ids:
            dropped.append((smi, f'all {len(ids)} conformers failed MMFF convergence'))
            continue

        scores = []
        for c in converged_ids:
            O3A = rdMolAlign.GetO3A(target, query_mol, prbCid=c)
            O3A.Align()
            scores.append(1.0 - rdShapeHelpers.ShapeTanimotoDist(
                target, query_mol, confId1=c,
            ))
        hits.append((target, max(scores)))

    if dropped:
        print(f'WARNING: {len(dropped)} library molecule(s) produced no usable '
              f'conformer and were dropped:')
        for smi, reason in dropped:
            print(f'  {smi}: {reason}')
    return sorted(hits, key=lambda x: x[1], reverse=True)


def ecfp_tanimoto(mol1, mol2):
    # Radius 2 / 2048 bits: compare only fingerprints built with identical settings
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    return DataStructs.TanimotoSimilarity(gen.GetFingerprint(mol1), gen.GetFingerprint(mol2))

# These thresholds are repository starting defaults only; calibrate both on a
# task-relevant active/decoy or retrieval benchmark before making decisions.
def scaffold_hop_candidates(query_mol, library, shape_threshold=0.7,
                            ecfp_threshold=0.5, n_conf=20, seed=42):
    shape_hits = shape_search_ensemble(query_mol, library, n_conf=n_conf, seed=seed)
    candidates = []
    for target, shape_score in shape_hits:
        if shape_score >= shape_threshold:
            ecfp_sim = ecfp_tanimoto(query_mol, target)
            if ecfp_sim < ecfp_threshold:
                candidates.append((target, shape_score, ecfp_sim))
    return candidates


def main():
    ap = argparse.ArgumentParser(description='Conformer-ensemble Open3DAlign shape search and scaffold-hop filter')
    ap.add_argument('--query', required=True, help='query SMILES (single connected molecule)')
    ap.add_argument('--library', required=True, help='.smi file, one SMILES per line')
    ap.add_argument('--n-conf', type=int, default=20)
    ap.add_argument('--seed', type=int, default=42)
    ap.add_argument('--scaffold-hop', action='store_true',
                    help='keep only shape >= --shape-threshold and ECFP4 < --ecfp-threshold')
    ap.add_argument('--shape-threshold', type=float, default=0.7)
    ap.add_argument('--ecfp-threshold', type=float, default=0.5)
    a = ap.parse_args()

    query = Chem.MolFromSmiles(a.query)
    if query is None or len(Chem.GetMolFrags(query)) > 1:
        raise SystemExit('query must be one valid connected SMILES')
    query = Chem.AddHs(query)
    params = AllChem.ETKDGv3()
    params.randomSeed = a.seed
    if AllChem.EmbedMolecule(query, params) != 0:
        raise SystemExit('query embedding failed')

    library = []
    for line in open(a.library, encoding='utf-8'):
        if line.strip():
            mol = Chem.MolFromSmiles(line.split()[0])
            if mol is None:
                print(f'WARNING: unparsable SMILES skipped: {line.strip()}')
            else:
                library.append(mol)

    if a.scaffold_hop:
        for target, shape, ecfp in scaffold_hop_candidates(
                query, library, a.shape_threshold, a.ecfp_threshold, a.n_conf, a.seed):
            print(f'{Chem.MolToSmiles(Chem.RemoveHs(target))}	shape={shape:.3f}	ecfp4={ecfp:.3f}')
    else:
        for target, shape in shape_search_ensemble(query, library, n_conf=a.n_conf, seed=a.seed):
            print(f'{Chem.MolToSmiles(Chem.RemoveHs(target))}	shape={shape:.3f}')


if __name__ == '__main__':
    main()
