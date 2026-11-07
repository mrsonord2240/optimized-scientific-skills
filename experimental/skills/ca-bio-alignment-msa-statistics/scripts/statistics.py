#!/usr/bin/env python3
"""Write selected source alignment statistics to JSON."""
import argparse
import json
import math
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'examples'))
from msa_utils import load_alignment, pick_background, check_alphabet, is_nucleotide
from identity_matrix import identity_matrix_vectorized
from conservation_profile import column_conservation, average_conservation
from entropy_analysis import shannon_entropy, information_content
from gap_statistics import gap_statistics
from substitution_counts import substitution_counts, transition_transversion
from alignment_scores import alignment_score

def clean(value):
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if isinstance(value, (list, tuple)):
        return [clean(x) for x in value]
    if isinstance(value, dict):
        return {k: clean(v) for k, v in value.items()}
    return value

def main():
    p = argparse.ArgumentParser()
    p.add_argument('operation', choices=['identity', 'columns', 'scores'])
    p.add_argument('input')
    p.add_argument('output')
    p.add_argument('--method', choices=['pid1', 'pid2', 'pid3', 'pid4'], default='pid4')
    p.add_argument('--min-occupancy', type=float, default=0.5)
    args = p.parse_args()
    if not 0 <= args.min_occupancy <= 1:
        p.error('min occupancy must be between 0 and 1')
    a = load_alignment(args.input)
    result = {'ids': [r.id for r in a], 'columns': a.get_alignment_length()}
    if args.operation == 'identity':
        result.update(method=args.method, matrix=identity_matrix_vectorized(a, args.method).tolist())
    elif args.operation == 'columns':
        background, label = pick_background(a)
        result.update(background=label, excluded=check_alphabet(a, background),
                      conservation=[column_conservation(a, i, min_occupancy=args.min_occupancy) for i in range(a.get_alignment_length())],
                      average_conservation=average_conservation(a, min_occupancy=args.min_occupancy),
                      entropy=[shannon_entropy(a[:, i]) for i in range(a.get_alignment_length())],
                      information=[information_content(a[:, i], background) for i in range(a.get_alignment_length())],
                      gaps=gap_statistics(a))
    else:
        counts = substitution_counts(a)
        result.update(substitutions={'/'.join(pair): n for pair, n in counts.items()}, flat_gap_score=alignment_score(a))
        if is_nucleotide(a):
            result['transitions_transversions_other'] = transition_transversion(counts)
    Path(args.output).write_text(json.dumps(clean(result), indent=2, allow_nan=False) + '\n', encoding='utf-8')

if __name__ == '__main__':
    main()
