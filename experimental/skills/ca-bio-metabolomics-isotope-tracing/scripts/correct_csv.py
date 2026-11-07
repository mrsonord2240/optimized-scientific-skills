"""Correct one ordered 13C intensity vector with the source IsoCor helper."""
import argparse
import csv
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'examples'))
from isotope_correction import correct_isocor

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('input_csv')
parser.add_argument('formula')
parser.add_argument('purity', type=float)
parser.add_argument('output_csv')
args = parser.parse_args()
if not 0 < args.purity <= 1:
    raise ValueError('purity must be in (0, 1]')
with open(args.input_csv, newline='', encoding='utf-8') as handle:
    rows = list(csv.DictReader(handle))
raw = [float(row['intensity']) for row in rows]
if len(raw) < 2 or not all(math.isfinite(x) and x >= 0 for x in raw) or sum(raw) <= 0:
    raise ValueError('intensities must be finite, nonnegative and nonzero in total')
areas, fractions, enrichment = correct_isocor(raw, args.formula, [1 - args.purity, args.purity])
with open(args.output_csv, 'w', newline='', encoding='utf-8') as handle:
    writer = csv.writer(handle)
    writer.writerow(['isotopologue', 'corrected_area', 'fraction', 'mean_enrichment'])
    for index, (area, fraction) in enumerate(zip(areas, fractions, strict=True)):
        writer.writerow([f'M+{index}', area, fraction, enrichment])
