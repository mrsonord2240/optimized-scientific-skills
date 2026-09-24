'''Shared helpers for the msa-statistics examples: gap/case normalisation, alphabet check, backgrounds.

MAFFT writes nucleotide alignments in lower case, hmmalign/Stockholm/A2M write '.' gaps and lower-case
insert columns. Every statistic below compares against '-' and counts letters case-sensitively, so
un-normalised input gives silently wrong numbers (checked: Ti/Tv 0.00 instead of 1.21, DNA information
content 29.9 bits instead of 2.0). Normalise first, then run check_alphabet().
'''
# Reference: biopython 1.83+ (checked on 1.88), numpy 1.26+ | Verify API if version differs

import os
import sys
from collections import Counter

from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')
FORMAT_BY_EXTENSION = {
    '.sto': 'stockholm', '.stk': 'stockholm', '.stockholm': 'stockholm',
    '.aln': 'clustal', '.clw': 'clustal', '.clustal': 'clustal',
    '.phy': 'phylip-relaxed', '.phylip': 'phylip-relaxed',
    '.nex': 'nexus', '.nexus': 'nexus',
}

# Robinson & Robinson 1991 PNAS 88:8880 (values as tabulated by NCBI BLAST); sums to 1.0.
ROBINSON_BACKGROUND = {
    'A': 0.07805, 'R': 0.05129, 'N': 0.04487, 'D': 0.05364, 'C': 0.01925,
    'Q': 0.04264, 'E': 0.06295, 'G': 0.07377, 'H': 0.02199, 'I': 0.05142,
    'L': 0.09019, 'K': 0.05744, 'M': 0.02243, 'F': 0.03856, 'P': 0.05203,
    'S': 0.07120, 'T': 0.05841, 'W': 0.01330, 'Y': 0.03216, 'V': 0.06441,
}
DNA_UNIFORM = {'A': 0.25, 'C': 0.25, 'G': 0.25, 'T': 0.25}


def example_path(default_name):
    '''argv[1] when given, otherwise the tiny alignment shipped in examples/data.'''
    return sys.argv[1] if len(sys.argv) > 1 else os.path.join(DATA_DIR, default_name)


def guess_format(path):
    '''AlignIO format from the extension (see FORMAT_BY_EXTENSION), 'fasta' otherwise. Pass fmt= to
    load_alignment for anything else. hmmalign A2M is ragged and cannot be read by AlignIO (see alignment-io).'''
    return FORMAT_BY_EXTENSION.get(os.path.splitext(path)[1].lower(), 'fasta')


def is_nucleotide(alignment, min_fraction=0.9):
    '''True when >= min_fraction of the non-gap characters are A/C/G/T/U/N (any case), or when >= 50% are
    and the rest are IUPAC nucleotide ambiguity codes (R Y S W K M B D H V): DNA with 12% R/Y/S/W/K/M is
    still DNA, while protein has only ~30% A/C/G/T/N and ~66% of the full IUPAC set, so it never passes.
    Check the result with pick_background()'s label, and override by choosing the background yourself.'''
    text = ''.join(str(r.seq) for r in alignment).upper().replace('-', '').replace('.', '').replace('~', '')
    if not text:
        return False
    core = sum(text.count(c) for c in 'ACGTUN') / len(text)
    ambiguity = sum(text.count(c) for c in 'RYSWKMBDHV') / len(text)
    return core >= min_fraction or (core >= 0.5 and core + ambiguity >= 0.98)


def normalize_alignment(alignment, upper=True, u_to_t=False):
    '''Copy of `alignment` with '.' and '~' gaps turned into '-', letters upper-cased (upper=True)
    and, for RNA (u_to_t=True), U turned into T. upper=False keeps case (A2M/A3M: lower case marks insert
    columns; hmmalign A2M is ragged and cannot be loaded by AlignIO, see alignment/alignment-io).'''
    records = []
    for record in alignment:
        seq = str(record.seq).replace('.', '-').replace('~', '-')
        if upper:
            seq = seq.upper()
        if u_to_t:
            seq = seq.replace('U', 'T').replace('u', 't')
        records.append(SeqRecord(Seq(seq), id=record.id, name=record.name, description=record.description))
    return MultipleSeqAlignment(records)


def check_alphabet(alignment, alphabet, label='alignment'):
    '''Count (and warn about) non-gap letters outside `alphabet` (e.g. a background dict) in an
    already normalised alignment. Statistics that need a fixed alphabet drop these letters.'''
    letters = Counter(''.join(str(r.seq) for r in alignment))
    outside = {c: n for c, n in letters.items() if c != '-' and c not in alphabet}
    if outside:
        total = sum(n for c, n in letters.items() if c != '-')
        print(f'WARNING {label}: {sum(outside.values())} of {total} residues are outside the alphabet '
              f'and are ignored by alphabet-bound statistics: {dict(sorted(outside.items()))}', file=sys.stderr)
    return outside


def pick_background(alignment):
    '''(background, label): DNA_UNIFORM for a nucleotide alignment, ROBINSON_BACKGROUND otherwise.'''
    if is_nucleotide(alignment):
        return DNA_UNIFORM, 'DNA (uniform background)'
    return ROBINSON_BACKGROUND, 'protein (Robinson 1991 background)'


def load_alignment(path, fmt=None):
    '''Read and normalise an alignment; format is guessed from the extension when fmt is None.
    RNA (U) is converted to T when the alignment is nucleotide.'''
    aln = AlignIO.read(path, fmt or guess_format(path))
    return normalize_alignment(aln, u_to_t=is_nucleotide(aln))
