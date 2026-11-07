'''Pairwise structural alignment with TMalign / USalign and parse the structured output.

Returns both TM-scores (normalised by chain 1 and by chain 2), RMSD, alignment length, and
optionally the superposed mobile structure. Uses -outfmt 2 for tabular output that is
robust to TM-align release versions.
'''
# Reference: TM-align 20220412+, US-align 20231222+ (checked on TM-align 20240303, US-align 20241108)
# | Verify CLI flags if version differs

import subprocess
import sys

def tm_align(reference_pdb, mobile_pdb, output_prefix=None, multimer=False):
    '''output_prefix is a prefix, not a file name: the tools write <output_prefix>.pdb (full-atom
    mobile structure superposed onto the reference) plus PyMOL .pml scripts.'''
    binary = 'USalign' if multimer else 'TMalign'
    cmd = [binary, mobile_pdb, reference_pdb, '-outfmt', '2']
    if multimer:
        cmd += ['-mm', '1', '-ter', '0']
    if output_prefix:
        cmd += ['-o', output_prefix]
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    parsed = parse_outfmt2(result.stdout)
    if parsed is None:
        # TM-align exits 0 with no data row for a <3-residue chain or a file without protein atoms
        message = ' '.join(l for l in (result.stdout + result.stderr).splitlines()
                           if l.strip() and not l.startswith('#'))
        raise RuntimeError(f'{binary} produced no alignment for {mobile_pdb} vs {reference_pdb}: {message}')
    if output_prefix:
        parsed['superposed_pdb'] = f'{output_prefix}.pdb'
    return parsed

def parse_outfmt2(stdout):
    '''Parse TMalign / USalign -outfmt 2 tabular output.

    Header columns: PDBchain1 PDBchain2 TM1 TM2 RMSD ID1 ID2 IDali L1 L2 Lali
    '''
    for line in stdout.splitlines():
        if line.startswith('#') or not line.strip():
            continue
        fields = line.split()
        if len(fields) < 11:
            continue
        return {
            'pdb1': fields[0],
            'pdb2': fields[1],
            'tm1': float(fields[2]),
            'tm2': float(fields[3]),
            'rmsd': float(fields[4]),
            'id1': float(fields[5]),
            'id2': float(fields[6]),
            'idali': float(fields[7]),
            'length1': int(fields[8]),
            'length2': int(fields[9]),
            'length_align': int(fields[10]),
        }
    return None

def interpret_tmscore(tm, length=None):
    '''Fold-level reading of a TM-score; fold similarity is not proof of homology.
    Pass the shorter chain length to apply the short-chain guard.'''
    if length is not None and length < 60:
        return 'chain < 60 residues: TM-score is unreliable, do not call a fold'
    if tm > 0.8:
        return 'very similar topology'
    if tm > 0.5:
        return 'same fold'
    if tm > 0.2:
        return 'weak structural similarity'
    return 'statistically random'

if __name__ == '__main__':
    # usage: python tm_align_pairwise.py REFERENCE.pdb MOBILE.pdb
    result = tm_align(sys.argv[1], sys.argv[2], output_prefix='superposed')
    tm_whole = min(result['tm1'], result['tm2'])  # normalised by the longer chain: whole-chain fold call
    print(f'TM-score normalised by chain 1 (mobile): {result["tm1"]:.3f}')
    print(f'TM-score normalised by chain 2 (reference): {result["tm2"]:.3f}')
    print(f'RMSD: {result["rmsd"]:.2f} A over {result["length_align"]} aligned residues '
          f'(L1={result["length1"]}, L2={result["length2"]})')
    print(f'Sequence identity in alignment: {result["idali"]*100:.1f}%')
    print(f'Interpretation (smaller TM, {min(result["length1"], result["length2"])}-residue shorter chain): '
          f'{interpret_tmscore(tm_whole, min(result["length1"], result["length2"]))}')
    print(f'Superposed mobile structure: {result["superposed_pdb"]}')
