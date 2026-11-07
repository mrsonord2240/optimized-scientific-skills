#!/usr/bin/env python3
# Purpose: subprocess wrapper around the NCBI Datasets CLI: run `datasets summary` and parse the
#          JSON-lines into dicts; run `datasets download` and return the zip path.
#          Reference: NCBI Datasets CLI 18.37.0 (checked 2026-09-19) | Verify API if version differs
# Inputs:  --taxon NAME     taxon for `datasets summary genome taxon` (reference assemblies only)
#          --accession ACC  optional; also download this assembly (--include list, --out zip)
# Usage:   python scripts/datasets_wrapper.py --taxon "Escherichia coli" \
#              [--accession GCF_000005845.2 --out ecoli_k12.zip --include genome,gff3,protein]
#          or import: from datasets_wrapper import datasets_summary, datasets_download
import argparse
import subprocess
import json
from pathlib import Path


def datasets_summary(subcommand, *args):
    '''Run `datasets summary` and parse JSON-lines stdout.'''
    cmd = ['datasets', 'summary', subcommand, *args, '--as-json-lines']
    out = subprocess.run(cmd, capture_output=True, text=True, check=True)
    return [json.loads(line) for line in out.stdout.strip().split('\n') if line]


def datasets_download(subcommand, *args, out='dataset.zip', include=None):
    cmd = ['datasets', 'download', subcommand, *args, '--filename', out]
    if include:
        cmd += ['--include', ','.join(include)]
    subprocess.run(cmd, check=True)
    return Path(out)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--taxon', default='Escherichia coli')
    ap.add_argument('--accession', default=None)
    ap.add_argument('--out', default='dataset.zip')
    ap.add_argument('--include', default='genome,gff3,protein')
    a = ap.parse_args()

    genomes = datasets_summary('genome', 'taxon', a.taxon, '--reference')
    print(f'{len(genomes)} reference {a.taxon} assemblies')
    for g in genomes[:3]:
        acc = g.get('accession')
        n50 = g.get('assembly_stats', {}).get('contig_n50')  # snake_case JSON keys, not camelCase
        print(f'  {acc}  N50={n50}')

    if a.accession:
        datasets_download('genome', 'accession', a.accession,
                          out=a.out,
                          include=a.include.split(','))


if __name__ == '__main__':
    main()
