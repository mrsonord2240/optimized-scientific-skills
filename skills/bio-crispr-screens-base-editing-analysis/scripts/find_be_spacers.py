"""Find base-editor spacers that place a target base in the editor's window at a target codon.

Purpose: tile NGG-adjacent spacers for a BE variant screen; each spacer is annotated with target vs
bystander editable bases (reverse-strand spacers are converted back to forward-CDS coordinates).
Inputs:  CDS nucleotide sequence (string or file), target amino-acid number, target base (C=CBE, A=ABE),
         editor name (BE3, BE4max, eA3A-BE3, ABE7.10, ABE8.20, ABE8e, evoCDA-BE).
Output:  table: spacer, strand, spacer_start, target_positions, bystander_positions, n_bystanders.
Usage:   python find_be_spacers.py --cds cds.fa --target-aa 130 --target-base C --editor BE4max --out spacers.tsv
Import:  from find_be_spacers import find_be_spacers
Needs:   pandas, biopython.
"""
import sys
import pandas as pd
import re
from Bio.Seq import Seq

def find_be_spacers(cds_sequence, cds_protein_start, target_aa, target_base='C', editor='BE4max'):
    '''Find sgRNAs that place target_base in editor-specific window at target_aa.
    Returns spacers with bystander annotation.

    Args:
        cds_sequence: nucleotide CDS (translated frame 1)
        cds_protein_start: amino acid number of CDS start (usually 1)
        target_aa: amino acid number to install variant (e.g., 130 for residue 130)
        target_base: 'C' (CBE) or 'A' (ABE)
        editor: 'BE3', 'BE4max', 'eA3A-BE3', 'ABE7.10', 'ABE8.20', 'ABE8e', 'evoCDA-BE'

    Returns: DataFrame (sorted by n_bystanders) with spacer, strand, spacer_start,
             target_positions, bystander_positions, n_bystanders. Empty (same columns) if no
             PAM-adjacent window holds target_base. Input is upper-cased before the
             case-sensitive PAM search.
    '''
    # Editor-specific editing window (positions from PAM-distal end of spacer)
    window_by_editor = {
        'BE3': (4, 8),       'BE4max': (4, 8),    'eA3A-BE3': (5, 7),
        'ABE7.10': (4, 7),   'ABE8.20': (4, 8),   'ABE8e': (4, 8),     # SpABE8e matches CBE window (Richter 2020)
        'evoCDA-BE': (1, 9),
    }
    if editor not in window_by_editor:
        raise ValueError(f"editor={editor!r} not recognized; valid editors: {sorted(window_by_editor)}")
    window_lo, window_hi = window_by_editor[editor]
    cds_sequence = cds_sequence.upper()
    aa_index = target_aa - cds_protein_start  # 0-indexed in protein
    aa_start_nt = aa_index * 3                # nt offset in cds
    candidates = []
    spacer_len = 20
    pam_pattern = re.compile(r'(?=([ACGT]GG))')
    for strand, seq in [('+', cds_sequence), ('-', str(Seq(cds_sequence).reverse_complement()))]:
        for pam_match in pam_pattern.finditer(seq):
            pam_pos = pam_match.start()
            spacer_start = pam_pos - spacer_len
            if spacer_start < 0:
                continue
            spacer = seq[spacer_start:pam_pos]
            # Editor-specific window from PAM-distal end (1-indexed)
            # Find all editable bases in window
            edit_bases_in_window = []
            for i, b in enumerate(spacer[window_lo-1:window_hi], start=window_lo):
                if b == target_base:
                    edit_bases_in_window.append(i)
            if not edit_bases_in_window:
                continue
            # Annotate which edits hit the target_aa codon
            target_codon_start = aa_start_nt
            target_codon_end = target_codon_start + 3
            target_position_in_spacer = []
            for i in edit_bases_in_window:
                pos_in_seq = spacer_start + i - 1  # 0-indexed position within `seq` (strand-specific)
                # For strand '-', `seq` is the reverse complement of cds_sequence; convert
                # back to forward-CDS coordinates before comparing against target_codon_start/
                # end, which are always forward-strand. Without this, reverse-strand spacers
                # silently misattribute target vs bystander (verified: a hand-constructed
                # reverse-strand case with a known on-target C was called "bystander" by the
                # unconverted math, and the audit's own random-CDS run produced an on-target
                # call 65nt from the true codon).
                genomic_pos = (len(cds_sequence) - 1 - pos_in_seq) if strand == '-' else pos_in_seq
                if target_codon_start <= genomic_pos < target_codon_end:
                    target_position_in_spacer.append(i)
            bystander_positions = [i for i in edit_bases_in_window if i not in target_position_in_spacer]
            candidates.append({
                'spacer': spacer,
                'strand': strand,
                'spacer_start': spacer_start,
                'target_positions': target_position_in_spacer,
                'bystander_positions': bystander_positions,
                'n_bystanders': len(bystander_positions),
            })
    cols = ['spacer', 'strand', 'spacer_start', 'target_positions', 'bystander_positions', 'n_bystanders']
    if not candidates:
        return pd.DataFrame(columns=cols)
    return pd.DataFrame(candidates, columns=cols).sort_values('n_bystanders')


if __name__ == '__main__':
    import argparse
    import os
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--cds', required=True, help='CDS nucleotide sequence, or a path to a file holding it (FASTA header ignored)')
    ap.add_argument('--target-aa', type=int, required=True, help='amino acid number to install the variant at')
    ap.add_argument('--cds-protein-start', type=int, default=1)
    ap.add_argument('--target-base', default='C', choices=['C', 'A'], help="'C' (CBE) or 'A' (ABE)")
    ap.add_argument('--editor', default='BE4max')
    ap.add_argument('--out', help='write TSV here (default: print)')
    a = ap.parse_args()
    cds = a.cds
    if os.path.exists(cds):
        cds = ''.join(l.strip() for l in open(cds) if not l.startswith('>'))
    res = find_be_spacers(cds, a.cds_protein_start, a.target_aa, a.target_base, a.editor)
    if a.out:
        res.to_csv(a.out, sep='\t', index=False)
    else:
        print(res.to_string(index=False))
    n0 = int((res.n_bystanders == 0).sum()) if len(res) else 0
    print(f'{len(res)} candidate spacers; {n0} with zero bystanders', file=sys.stderr)
