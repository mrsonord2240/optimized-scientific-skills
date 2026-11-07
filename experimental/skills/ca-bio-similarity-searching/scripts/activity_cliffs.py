from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator

def activity_cliffs(df, sim_threshold=0.85, activity_gap=2.0, activity_col='pIC50',
                    include_chirality=True):
    generator = rdFingerprintGenerator.GetMorganGenerator(
        radius=2, fpSize=2048, includeChirality=include_chirality)
    mols = [Chem.MolFromSmiles(s) for s in df['smiles']]
    if any(mol is None for mol in mols):
        raise ValueError('activity-cliff input contains invalid SMILES')
    fps = [generator.GetFingerprint(mol) for mol in mols]
    cliffs = []
    for i in range(len(fps)):
        sims = DataStructs.BulkTanimotoSimilarity(fps[i], fps[i+1:])
        for j_off, sim in enumerate(sims):
            j = i + 1 + j_off
            if sim >= sim_threshold:
                gap = abs(df[activity_col].iloc[i] - df[activity_col].iloc[j])
                if gap >= activity_gap:
                    cliffs.append((i, j, sim, gap))
    return cliffs
