#!/usr/bin/env python3
"""Inspect, convert, or slice rectangular alignments using Bio.AlignIO."""
import argparse
import json
from pathlib import Path
from Bio import AlignIO


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('operation', choices=['inspect', 'convert', 'slice'])
    parser.add_argument('input')
    parser.add_argument('input_format')
    parser.add_argument('output')
    parser.add_argument('output_format', nargs='?')
    parser.add_argument('--molecule-type', choices=['DNA', 'RNA', 'protein'])
    parser.add_argument('--start', type=int, default=0)
    parser.add_argument('--end', type=int)
    parser.add_argument('--row-start', type=int, default=0)
    parser.add_argument('--row-end', type=int)
    args = parser.parse_args()
    if Path(args.input).resolve() == Path(args.output).resolve():
        parser.error('input and output must differ')
    alignments = list(AlignIO.parse(args.input, args.input_format))
    if not alignments:
        parser.error('input contains no alignments')
    if args.operation == 'inspect':
        summary = [{'sequences': len(a), 'columns': a.get_alignment_length(),
                    'ids': [record.id for record in a]} for a in alignments]
        Path(args.output).write_text(json.dumps(summary, indent=2) + '\n', encoding='utf-8')
        return
    if not args.output_format:
        parser.error('output format is required')
    if args.output_format == 'nexus' and not args.molecule_type:
        parser.error('NEXUS requires an explicit molecule type')
    if args.operation == 'slice':
        if len(alignments) != 1:
            parser.error('slice requires exactly one alignment')
        alignment = alignments[0]
        end = args.end if args.end is not None else alignment.get_alignment_length()
        row_end = args.row_end if args.row_end is not None else len(alignment)
        if not (0 <= args.start < end <= alignment.get_alignment_length()
                and 0 <= args.row_start < row_end <= len(alignment)):
            parser.error('slice is empty or outside the alignment')
        alignments = [alignment[args.row_start:row_end, args.start:end]]
    if args.molecule_type:
        for alignment in alignments:
            for record in alignment:
                record.annotations['molecule_type'] = args.molecule_type
    AlignIO.write(alignments, args.output, args.output_format)


if __name__ == '__main__':
    main()
