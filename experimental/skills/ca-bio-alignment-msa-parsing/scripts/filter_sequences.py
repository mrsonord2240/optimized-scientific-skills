'''Filter the sequences (rows) of an alignment: by ID pattern, gap content, and/or duplicates.

Inputs : alignment file (fasta, or stockholm for .sto/.stk/.stockholm; override with --format)
Output : the kept rows, written in the same format (the ORIGINAL records, unchanged)
Usage  : python scripts/filter_sequences.py in.fasta out.fasta [--id-pattern REGEX]
                 [--max-gap-fraction 0.1] [--dedup]
Import : sys.path.insert(0, 'scripts'); from filter_sequences import filter_by_id, filter_by_gap_content, remove_duplicates

A filter that would remove every sequence raises ValueError instead of returning an empty alignment.
Gap fractions and duplicates are computed on the normalised alignment ('.' -> '-', upper case), so
AC-GT, AC.GT and ac-gt count as one sequence; column and record annotations of the input are kept.
'''
# Reference: biopython 1.83+ (checked on 1.88) | Verify API if version differs

import argparse
import os
import re
import sys

from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'examples'))
from msa_utils import guess_format, normalize_alignment  # noqa: E402


def _require_kept(kept, alignment, what):
    if not kept:
        raise ValueError(f'{what} removes all {len(alignment)} sequences')
    return MultipleSeqAlignment(kept, annotations=alignment.annotations,
                                column_annotations=alignment.column_annotations)


def filter_by_id(alignment, pattern):
    regex = re.compile(pattern)
    return _require_kept([r for r in alignment if regex.search(r.id)], alignment, f'pattern {pattern!r}')


def filter_by_gap_content(alignment, max_gap_fraction=0.1):
    fractions = [str(r.seq).count('-') / len(r.seq) for r in normalize_alignment(alignment)]
    kept = [r for r, f in zip(alignment, fractions) if f <= max_gap_fraction]
    return _require_kept(kept, alignment, f'max_gap_fraction={max_gap_fraction} (lowest fraction {min(fractions):.2f})')


def remove_duplicates(alignment):
    # Rows are compared after normalisation (AC-GT, AC.GT and ac-gt are one sequence); the original records are kept.
    seen, kept = set(), []
    for record, normalized in zip(alignment, normalize_alignment(alignment)):
        key = str(normalized.seq)
        if key not in seen:
            seen.add(key)
            kept.append(record)
    return _require_kept(kept, alignment, 'remove_duplicates')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('input')
    parser.add_argument('output')
    parser.add_argument('--format', help='alignment format (default: guessed from the extension)')
    parser.add_argument('--id-pattern', help='keep rows whose id matches this regex (re.search)')
    parser.add_argument('--max-gap-fraction', type=float, help='drop rows with a larger gap fraction')
    parser.add_argument('--dedup', action='store_true', help='drop rows identical after normalisation')
    args = parser.parse_args()

    fmt = args.format or guess_format(args.input)
    alignment = AlignIO.read(args.input, fmt)
    print(f'{args.input}: {len(alignment)} sequences x {alignment.get_alignment_length()} columns')
    if args.id_pattern:
        alignment = filter_by_id(alignment, args.id_pattern)
        print(f'  --id-pattern {args.id_pattern!r}: {len(alignment)} left')
    if args.max_gap_fraction is not None:
        alignment = filter_by_gap_content(alignment, args.max_gap_fraction)
        print(f'  --max-gap-fraction {args.max_gap_fraction}: {len(alignment)} left')
    if args.dedup:
        alignment = remove_duplicates(alignment)
        print(f'  --dedup: {len(alignment)} left')
    AlignIO.write(alignment, args.output, fmt)
    print(f'Wrote {len(alignment)} sequences to {args.output}')
