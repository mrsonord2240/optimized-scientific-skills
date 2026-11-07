"""Split a source-format MS-DIAL alignment table into metadata and sample values."""
import argparse
import csv
from pathlib import Path

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('input_tsv')
parser.add_argument('output_prefix')
args = parser.parse_args()
with open(args.input_tsv, newline='', encoding='utf-8-sig') as handle:
    rows = list(csv.reader(handle, delimiter='\t'))
if len(rows) < 5 or rows[0].count('Class') != 1:
    raise ValueError('expected four header rows and exactly one Class anchor')
anchor = rows[0].index('Class')
header = rows[4]
start = anchor + 1
if start >= len(header) or len(set(header)) != len(header) or 'Alignment ID' not in header:
    raise ValueError('invalid sample boundary, duplicate headers or missing Alignment ID')
body = rows[5:]
if any(len(row) != len(header) for row in body):
    raise ValueError('ragged data rows')
for row in body:
    for value in row[start:]:
        if value.strip():
            float(value)
prefix = Path(args.output_prefix)
for suffix, columns in (('.features.tsv', list(range(start))),
                        ('.intensity.tsv', [header.index('Alignment ID')] + list(range(start, len(header))))):
    with open(str(prefix) + suffix, 'w', newline='', encoding='utf-8') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        writer.writerow([header[i] for i in columns])
        writer.writerows([row[i] for i in columns] for row in body)
