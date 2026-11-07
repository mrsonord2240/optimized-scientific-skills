#!/usr/bin/env python3
"""Apply source default gap-based column and row filters."""
from pathlib import Path
import sys
from Bio import AlignIO

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'examples'))
from clean_alignment import remove_gappy_columns, filter_by_gap_content
from msa_utils import guess_format, normalize_alignment

if __name__ == '__main__':
    if len(sys.argv) != 3:
        raise SystemExit('usage: clean_columns.py input output')
    alignment = AlignIO.read(sys.argv[1], guess_format(sys.argv[1]))
    filtered = filter_by_gap_content(remove_gappy_columns(alignment, threshold=0.5), max_gap_fraction=0.2)
    AlignIO.write(normalize_alignment(filtered, upper=False), sys.argv[2], guess_format(sys.argv[2]))
