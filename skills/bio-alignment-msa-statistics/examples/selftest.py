'''Check the msa-statistics functions against values worked out by hand on the shipped tiny alignments.

Usage: python selftest.py     (prints OK per check; raises AssertionError on the first mismatch)

example_protein.fasta   s1 AC-DEF--   s2 -CGDEFHK   s3 ACGN-FHR   (terminal gaps, internal gap, one substitution pair)
example_dna.fasta       d1 acgtacgtac  d2 acgtgcgtat  d3 .cctacgnac  (lower case, '.' gap, N)
'''
import math

import numpy as np
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from capra_singh_jsd import capra_singh_score
from conservation_profile import column_conservation, conservation_profile
from entropy_analysis import information_content, shannon_entropy
from identity_matrix import average_identity, identity_matrix_vectorized, pairwise_identity
from kimura_protein_distance import kimura_protein_distance
from msa_utils import DNA_UNIFORM, ROBINSON_BACKGROUND, example_path, load_alignment, normalize_alignment
from pssm import pssm_with_pseudocounts
from substitution_counts import substitution_counts, transition_transversion


def close(a, b, tol=1e-4):
    return abs(a - b) <= tol


def ok(name):
    print(f'OK  {name}')


protein = load_alignment(example_path('example_protein.fasta'))
rows = [str(r.seq) for r in protein]

# Percent identity s1 vs s2: identical C,D,E,F = 4; aligned pairs 4; internal gap columns 1; ungapped lengths 5 and 7
for method, expected in [('pid1', 4 / 5), ('pid2', 4 / 4), ('pid3', 4 / 5), ('pid4', 4 / 6)]:
    assert close(pairwise_identity(rows[0], rows[1], method), expected), method
    assert close(identity_matrix_vectorized(protein, method)[0, 1], expected), method
ok('PID1-PID4 (function and vectorized) = 80.0 / 100.0 / 80.0 / 66.7 %')

assert math.isnan(pairwise_identity('AC--', '----', 'pid4')) and math.isnan(pairwise_identity('AC--', '--GT', 'pid2'))
matrix = identity_matrix_vectorized(MultipleSeqAlignment([SeqRecord(Seq('AC-D'), id='a'), SeqRecord(Seq('----'), id='b')]))
assert matrix[0, 0] == 1.0 and np.isnan(matrix[[0, 1, 1], [1, 0, 1]]).all()
m4 = identity_matrix_vectorized(protein, 'pid4')
avg, undefined = average_identity(m4)
assert close(avg, np.mean([m4[0, 1], m4[0, 2], m4[1, 2]])) and undefined == 0
ok('undefined identity (all-gap sequence) is NaN, average ignores it')

# Conservation per column [1, 1, 1, 2/3, 1, 1, 1, 1/2]; mean 0.8958 when every occupied column counts
cons = [column_conservation(protein, i, min_occupancy=0.0) for i in range(8)]
assert all(close(a, b) for a, b in zip(cons, [1, 1, 1, 2 / 3, 1, 1, 1, 1 / 2]))
assert close(np.mean(cons), 0.8958, 1e-3)
assert math.isnan(column_conservation(MultipleSeqAlignment([SeqRecord(Seq('A-'), id='a'), SeqRecord(Seq('--'), id='b'),
                                                               SeqRecord(Seq('--'), id='c')]), 0, min_occupancy=0.5))
assert math.isnan(column_conservation(protein, 0, min_occupancy=0.9))  # column 0 has 2 of 3 residues
prof = conservation_profile(protein, window=2, min_occupancy=0.0)  # centred window i-1..i+1
assert close(prof[1], np.mean(cons[0:3])) and close(prof[7], np.mean(cons[6:8]))
ok('conservation, min_occupancy NaN rule, centred profile window')

# Shannon entropy (bits): only column 3 (D,D,N: 0.9183) and column 7 (-,K,R: 1.0) are variable
ent = [shannon_entropy(protein[:, i]) for i in range(8)]
assert all(close(a, b) for a, b in zip(ent, [0, 0, 0, 0.9183, 0, 0, 0, 1.0]))
ok('Shannon entropy')

# Information content: KL against a background; letters outside it are dropped, not scored at 1e-9
assert close(information_content('AAAA', DNA_UNIFORM), 2.0)
assert close(information_content('AAAAN', DNA_UNIFORM), 2.0) and close(information_content('AAAAX-', ROBINSON_BACKGROUND),
                                                                       math.log2(1 / ROBINSON_BACKGROUND['A']))
assert information_content('NNNN', DNA_UNIFORM) == 0.0
ok('information content: unknown letters dropped and renormalised, DNA IC of a conserved column = 2.0 bits')

# Pseudocount PSSM column 1 (C,C,C): log2((3 + bg_C) / (3 + 1) / bg_C); ambiguity letters do not enter n
pssm = pssm_with_pseudocounts(protein, ROBINSON_BACKGROUND)
bg = ROBINSON_BACKGROUND['C']
assert close(pssm[1]['C'], math.log2((3 + bg) / 4 / bg))
ambiguous = normalize_alignment(MultipleSeqAlignment([SeqRecord(Seq('CX'), id='a'), SeqRecord(Seq('CC'), id='b')]))
assert close(pssm_with_pseudocounts(ambiguous, ROBINSON_BACKGROUND)[1]['C'], math.log2((1 + bg) / 2 / bg))
ok('PSSM with background-weighted pseudocount')

# Kimura distance s2 vs s3: 6 aligned pairs, 4 identical, p = 1/3, d = -ln(1 - 1/3 - 0.2/9) = 0.4393
assert close(kimura_protein_distance(rows[1], rows[2]), 0.4393)
ok('Kimura distance 0.4393')

# Substitution counts (protein): D/N in column 3 for two pairs, K/R once
subs = substitution_counts(protein)
assert subs == {('D', 'N'): 2, ('K', 'R'): 1}, subs
ok('substitution counts D-N 2, K-R 1')

# Nucleotide: lower case + '.' gap + N are normalised; 4 transitions (col 4 a/g x2, col 9 c/t x2), 2 transversions (col 2 g/c x2)
dna = load_alignment(example_path('example_dna.fasta'))
assert [str(r.seq) for r in dna][2] == '-CCTACGNAC'
ti, tv, other = transition_transversion(substitution_counts(dna))
assert (ti, tv, other) == (4, 2, 2), (ti, tv, other)
assert close(information_content(dna[:, 4], DNA_UNIFORM), 2 - 0.9183) and close(information_content(dna[:, 7], DNA_UNIFORM), 2.0)
ok('DNA: Ti/Tv 4/2 (N pairs excluded), IC 1.082 and 2.000 bits')

# Case, '.' and '~' gaps and RNA U must not change any statistic
messy = MultipleSeqAlignment([SeqRecord(Seq(s.lower().replace('-', '.')), id=r.id) for s, r in zip(rows, protein)])
messy = normalize_alignment(messy)
assert [str(r.seq) for r in messy] == rows
rna = normalize_alignment(MultipleSeqAlignment([SeqRecord(Seq('acgu~'), id='r')]), u_to_t=True)
assert str(rna[0].seq) == 'ACGT-'
ok('normalisation: lower case, "." / "~" gaps, U -> T')

# Sum-of-pairs conventions are defined in SKILL.md; here only the capra_singh_score shape/NaN guard
scores = capra_singh_score(protein)
assert len(scores) == 8 and all(np.isfinite(scores))
ok('capra_singh_score returns one finite score per column')
print('\nAll checks passed.')
