#!/usr/bin/env python3
"""Align exactly two FASTA records and save the first optimal alignment."""
import argparse
import json
from pathlib import Path
from Bio import SeqIO
from Bio.Align import PairwiseAligner, substitution_matrices
from Bio.Seq import Seq

def main():
    p = argparse.ArgumentParser()
    p.add_argument('input')
    p.add_argument('output')
    p.add_argument('--mode', choices=['global', 'local', 'semiglobal'], default='global')
    p.add_argument('--type', choices=['dna', 'protein'], default='dna')
    p.add_argument('--reverse-complement', action='store_true')
    a = p.parse_args()
    records = list(SeqIO.parse(a.input, 'fasta'))
    if len(records) != 2:
        p.error('input must contain exactly two FASTA records')
    seqs = [str(r.seq).strip().upper() for r in records]
    if not all(seqs):
        p.error('sequences must be non-empty')
    if a.reverse_complement and a.type != 'dna':
        p.error('reverse complement requires DNA')
    if a.type == 'dna':
        seqs = [s.replace('U', 'T') for s in seqs]
        aligner = PairwiseAligner(mode='local' if a.mode == 'local' else 'global', match_score=2, mismatch_score=-1, open_gap_score=-10, extend_gap_score=-0.5)
    else:
        aligner = PairwiseAligner(mode='local' if a.mode == 'local' else 'global', substitution_matrix=substitution_matrices.load('BLOSUM62'), open_gap_score=-11, extend_gap_score=-1)
    if a.mode == 'semiglobal':
        aligner.end_gap_score = 0.0
    strand = 'forward'
    if a.reverse_complement:
        reverse = str(Seq(seqs[1]).reverse_complement())
        if aligner.score(seqs[0], reverse) > aligner.score(*seqs):
            seqs[1] = reverse
            strand = 'reverse-complement'
    alignments = aligner.align(*seqs)
    alignment = alignments[0]
    counts = alignment.counts()
    denom = counts.identities + counts.mismatches
    result = {'ids': [r.id for r in records], 'mode': a.mode, 'type': a.type, 'query_strand': strand,
              'score': alignment.score, 'alignment': str(alignment), 'coordinates': alignment.coordinates.tolist(),
              'identities': counts.identities, 'mismatches': counts.mismatches, 'gaps': counts.gaps,
              'pid2': counts.identities / denom if denom else None, 'scoring': str(aligner)}
    Path(a.output).write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')

if __name__ == '__main__':
    main()
