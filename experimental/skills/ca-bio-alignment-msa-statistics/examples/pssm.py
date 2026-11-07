'''Position-specific scoring matrix (PSSM) with background-weighted pseudocounts.

Reference: BioPython 1.83+ (checked on 1.88) | Verify API if version differs
Reference: Henikoff JG & Henikoff S 1996 Bioinf 12:135-143 (data-dependent pseudocounts).

Without pseudocounts, log-odds against background diverge to negative infinity at any column missing
a residue. Here `pseudocount` is the total pseudocount per column (default 1.0), spread over residues
in proportion to the background; HMMER uses Dirichlet mixtures for finer control.

Input must be normalised (upper case, '-' gaps); letters outside the background (X, B, Z, N, ...) are
dropped from the counts and from the column total. The default background is chosen by alphabet
(protein: Robinson & Robinson 1991; nucleotide: uniform), because A/C/G/T are in both alphabets and a
protein background would otherwise run silently on DNA.
'''
import math
from collections import Counter

from msa_utils import check_alphabet, example_path, load_alignment, pick_background


def pssm_with_pseudocounts(alignment, background=None, pseudocount=1.0):
    if background is None:
        background = pick_background(alignment)[0]
    pssm = []
    for col_idx in range(alignment.get_alignment_length()):
        column = [c for c in alignment[:, col_idx] if c in background]
        n = len(column)
        counts = Counter(column)
        pssm.append({
            residue: math.log2(((counts.get(residue, 0) + pseudocount * background[residue]) / (n + pseudocount))
                                / background[residue])
            for residue in background
        })
    return pssm


def score_site(sequence, pssm):
    return sum(pssm[i].get(residue, 0) for i, residue in enumerate(sequence))


if __name__ == '__main__':
    alignment = load_alignment(example_path('example_protein.fasta'))
    background, label = pick_background(alignment)
    print(f'Treating as {label}')
    check_alphabet(alignment, background, 'PSSM')
    pssm = pssm_with_pseudocounts(alignment, background)
    print(f'PSSM with {len(pssm)} positions')
    for i, col_scores in enumerate(pssm[:5]):
        top = sorted(col_scores.items(), key=lambda kv: -kv[1])[:3]
        print(f'  Position {i}: {", ".join(f"{r}={s:+.2f}" for r, s in top)}')
