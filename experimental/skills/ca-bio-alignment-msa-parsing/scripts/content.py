#!/usr/bin/env python3
"""Summarize an existing alignment without executing external tools."""
import argparse
from collections import Counter
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'examples'))
from msa_utils import load_alignment
from consensus_sequence import consensus_sequence
from find_conserved import find_conserved_positions

def main():
    p = argparse.ArgumentParser()
    p.add_argument('input')
    p.add_argument('output')
    p.add_argument('--consensus-threshold', type=float, default=0.5)
    p.add_argument('--conserved-threshold', type=float, default=0.8)
    args = p.parse_args()
    if not all(0 <= x <= 1 for x in (args.consensus_threshold, args.conserved_threshold)):
        p.error('thresholds must be in [0, 1]')
    a = load_alignment(args.input)
    result = {'ids': [r.id for r in a], 'columns': a.get_alignment_length(),
              'gaps_per_sequence': {r.id: str(r.seq).count('-') for r in a},
              'column_counts': [dict(Counter(a[:, i])) for i in range(a.get_alignment_length())],
              'consensus': consensus_sequence(a, threshold=args.consensus_threshold),
              'conserved': find_conserved_positions(a, threshold=args.conserved_threshold)}
    Path(args.output).write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')

if __name__ == '__main__':
    main()
