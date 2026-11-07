#!/usr/bin/env python3
"""
Grade ATAC-seq QC metrics against ENCODE-style thresholds; write one row per sample.
Reference: pandas 2.2+

Usage:
    python aggregate_qc.py metrics.json --output atac_qc_mqc.tsv [--sample NAME]

Input JSON is either {metric: value} (one sample, named by --sample, default "sample") or
{sample: {metric: value}}. Extra keys are ignored with a stderr note. null/NaN/absent -> "NA".
A non-numeric value stops the run with exit 2.

Output: sample-wide TSV with MultiQC custom-content headers (name it *_mqc.tsv): the seven
metric values, a <metric>_grade column each, and overall_grade (FAIL > INCOMPLETE > WARN > PASS).
"""

import argparse
import json
import math
import sys

# metric: (acceptable bound, ideal bound); mt_fraction is inverted (lower is better)
ENCODE_THRESHOLDS = {
    'nuclear_reads_M': (25, 50),
    'mt_fraction': (0.5, 0.05),
    'NRF': (0.7, 0.9),
    'PBC1': (0.7, 0.9),
    'PBC2': (1.0, 3.0),
    'TSS_enrichment': (5.0, 7.0),
    'FRiP': (0.2, 0.3),
}
INVERTED = {'mt_fraction'}


def to_number(metric, value):
    if value is None:
        return None
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f'metric {metric!r} must be a number or null, got {value!r}')
    return None if math.isnan(value) else float(value)


def grade(value, acceptable, ideal, inverted=False):
    if value is None:
        return 'NA'
    if inverted:
        return 'FAIL' if value > acceptable else ('PASS' if value <= ideal else 'WARN')
    return 'FAIL' if value < acceptable else ('PASS' if value >= ideal else 'WARN')


def overall(grades):
    if 'FAIL' in grades:
        return 'FAIL'
    if 'NA' in grades:
        return 'INCOMPLETE'
    return 'WARN' if 'WARN' in grades else 'PASS'


def report(samples, out_tsv=None):
    """samples: {sample: {metric: value}} -> list of row dicts (also written as TSV)."""
    rows = []
    for name, metrics in samples.items():
        extra = sorted(set(metrics) - set(ENCODE_THRESHOLDS))
        if extra:
            print(f'{name}: ignoring ungraded keys {extra}', file=sys.stderr)
        row, grades = {'sample': name}, []
        for m, (acc, ideal) in ENCODE_THRESHOLDS.items():
            v = to_number(m, metrics.get(m))
            g = grade(v, acc, ideal, m in INVERTED)
            row[m], row[f'{m}_grade'] = ('NA' if v is None else v), g
            grades.append(g)
        row['overall_grade'] = overall(grades)
        rows.append(row)
    cols = list(rows[0]) if rows else ['sample']
    header = ["# id: 'atac_encode_qc'", "# section_name: 'ATAC-seq ENCODE QC'", "# plot_type: 'table'"]
    lines = header + ['\t'.join(cols)] + ['\t'.join(str(r[c]) for c in cols) for r in rows]
    text = '\n'.join(lines) + '\n'
    if out_tsv:
        with open(out_tsv, 'w', encoding='utf-8', newline='\n') as f:
            f.write(text)
    else:
        sys.stdout.write(text)
    return rows


if __name__ == '__main__':
    p = argparse.ArgumentParser(description='Grade ATAC-seq QC metrics against ENCODE-style thresholds')
    p.add_argument('metrics', help='JSON file (or "-" for stdin)')
    p.add_argument('--output', help='Output TSV (default: stdout)')
    p.add_argument('--sample', default='sample', help='Sample name for a flat {metric: value} input')
    a = p.parse_args()
    data = json.load(sys.stdin if a.metrics == '-' else open(a.metrics, encoding='utf-8'))
    if not isinstance(data, dict):
        sys.exit('metrics JSON must be an object')
    samples = data if data and all(isinstance(v, dict) for v in data.values()) else {a.sample: data}
    try:
        report(samples, a.output)
    except ValueError as e:
        print(f'error: {e}', file=sys.stderr)
        sys.exit(2)
