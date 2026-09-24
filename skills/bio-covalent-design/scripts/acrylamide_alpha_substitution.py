# Count non-warhead neighbors on the alpha carbon of the first acrylamide match.
# A structural feature only -- not a LUMO estimate or a reactivity prediction.
# Input: one or more SMILES as arguments.
# Usage: python scripts/acrylamide_alpha_substitution.py 'C=CC(=O)N1CCCCC1' 'C=C(C)C(=O)N1CCCCC1'
# Prints "<smiles>\t<count or None>" (None = unparsable SMILES or no acrylamide).
import sys

from rdkit import Chem


def acrylamide_alpha_substitution_count(smi):
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None
    acryl_pat = Chem.MolFromSmarts(
        '[CX3:1](=[OX1:2])([NX3:3])[CX3:4]=[CX3:5]'
    )
    matches = mol.GetSubstructMatches(acryl_pat, uniquify=True)
    if not matches:
        return None
    alpha_query_idx = next(
        atom.GetIdx() for atom in acryl_pat.GetAtoms()
        if atom.GetAtomMapNum() == 4
    )
    alpha_c = mol.GetAtomWithIdx(matches[0][alpha_query_idx])
    n_subs = len([n for n in alpha_c.GetNeighbors() if n.GetIdx() not in matches[0]])
    return n_subs


if __name__ == '__main__':
    for smiles in sys.argv[1:]:
        print(f'{smiles}\t{acrylamide_alpha_substitution_count(smiles)}')
