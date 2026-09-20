'''Pure-Python pairwise structure superposition with Bio.PDB.Superimposer.

Use when residue correspondence is known (e.g. apo vs holo of the same protein,
mutant vs wild-type). For unknown correspondence, prefer TMalign / USalign.

CA atoms are paired by (chain, residue number, insertion code), never by list position, and
only standard residues of the first model are used (residue.id[0] == ' '). That excludes
waters, ligands and metal ions: a Ca2+ ion's atom is also named 'CA' and a name-only filter
pairs it with a real residue. The script refuses when the structures are not co-numbered.
'''
# Reference: biopython 1.83+ (checked on 1.88) | Verify API if version differs

from Bio.PDB import PDBParser, Superimposer, PDBIO


def ca_by_residue(structure):
    '''{(chain, resseq, icode): (resname, CA atom)} for standard residues of the first model.'''
    cas = {}
    for chain in structure[0]:
        for residue in chain:
            if residue.id[0] == ' ' and 'CA' in residue:
                cas[(chain.id, residue.id[1], residue.id[2])] = (residue.get_resname(), residue['CA'])
    return cas


def superpose_ca(reference_structure, mobile_structure, min_paired=0.5, max_mismatch=0.2):
    '''Superpose mobile onto reference in place; return (Superimposer, n_pairs, n_mismatched_names).

    Raises ValueError if fewer than min_paired of the shorter chain's residues share a
    (chain, number, icode) key, or if more than max_mismatch of the paired residue names differ
    (an offset numbering; point mutants stay far below this).
    '''
    reference = ca_by_residue(reference_structure)
    mobile = ca_by_residue(mobile_structure)
    keys = sorted(set(reference) & set(mobile))
    if not keys:
        raise ValueError('No residues share (chain, number, insertion code); the structures are not '
                         'co-numbered. Use TMalign / USalign, which find the correspondence themselves.')
    shorter = min(len(reference), len(mobile))
    if len(keys) < min_paired * shorter:
        raise ValueError(f'Only {len(keys)} of {shorter} CA atoms share a residue key; '
                         'use TMalign / USalign for structures that are not co-numbered.')
    n_mismatch = sum(reference[k][0] != mobile[k][0] for k in keys)
    if n_mismatch > max_mismatch * len(keys):
        raise ValueError(f'{n_mismatch} of {len(keys)} paired residues have different names; '
                         'the numbering is offset. Use TMalign / USalign.')
    sup = Superimposer()
    sup.set_atoms([reference[k][1] for k in keys], [mobile[k][1] for k in keys])
    sup.apply(list(mobile_structure.get_atoms()))
    return sup, len(keys), n_mismatch


if __name__ == '__main__':
    parser = PDBParser(QUIET=True)
    mobile_structure = parser.get_structure('mobile', 'mobile.pdb')
    reference_structure = parser.get_structure('reference', 'reference.pdb')

    sup, n_pairs, n_mismatch = superpose_ca(reference_structure, mobile_structure)

    print(f'Superposed {n_pairs} CA pairs ({n_mismatch} with different residue names)')
    print(f'RMSD: {sup.rms:.3f} A')
    print(f'Rotation matrix:\n{sup.rotran[0]}')
    print(f'Translation vector: {sup.rotran[1]}')

    io = PDBIO()
    io.set_structure(mobile_structure)
    io.save('mobile_superposed.pdb')
    print('\nSuperposed structure written to mobile_superposed.pdb')
