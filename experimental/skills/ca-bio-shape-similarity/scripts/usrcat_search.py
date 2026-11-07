import argparse
import csv
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

def descriptor(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None or mol.GetNumAtoms() == 0 or len(Chem.GetMolFrags(mol)) != 1:
        raise ValueError('SMILES must be one valid connected molecule')
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    if AllChem.EmbedMolecule(mol, params) != 0:
        raise ValueError('embedding failed')
    return rdMolDescriptors.GetUSRCAT(mol)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('query_smi')
    ap.add_argument('library_smi')
    ap.add_argument('output_csv')
    args = ap.parse_args()
    query = descriptor(open(args.query_smi, encoding='utf-8').readline().strip().split()[0])
    rows = []
    with open(args.library_smi, encoding='utf-8') as handle:
        for line_no, line in enumerate(handle, 1):
            if not line.strip():
                continue
            smi = line.split()[0]
            try:
                score = rdMolDescriptors.GetUSRScore(query, descriptor(smi))
                rows.append((line_no, smi, score, 'ok'))
            except (ValueError, RuntimeError) as exc:
                rows.append((line_no, smi, '', str(exc)))
    with open(args.output_csv, 'w', newline='', encoding='utf-8') as handle:
        writer = csv.writer(handle)
        writer.writerow(['line', 'smiles', 'usrcat_similarity', 'status'])
        writer.writerows(rows)

if __name__ == '__main__':
    main()
