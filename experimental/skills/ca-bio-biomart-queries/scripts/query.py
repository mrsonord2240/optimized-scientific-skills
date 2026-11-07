#!/usr/bin/env python3
"""Query one BioMart dataset and write a TSV using the source helper."""
import argparse
from pathlib import Path
from pybiomart import Server
from biomart_query import query_raw

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--dataset', required=True)
    p.add_argument('--attributes', required=True, help='comma-separated BioMart attribute names')
    p.add_argument('--filter', required=True)
    group = p.add_mutually_exclusive_group(required=True)
    group.add_argument('--value')
    group.add_argument('--values-file')
    p.add_argument('--output', required=True)
    p.add_argument('--host', default='http://www.ensembl.org')
    args = p.parse_args()
    attributes = [x.strip() for x in args.attributes.split(',') if x.strip()]
    values = ([x.strip() for x in Path(args.values_file).read_text(encoding='utf-8').splitlines() if x.strip()]
              if args.values_file else args.value)
    if not attributes or not values:
        p.error('attributes and filter values must not be empty')
    dataset = Server(host=args.host)['ENSEMBL_MART_ENSEMBL'][args.dataset]
    frame = query_raw(dataset, attributes=attributes, filters={args.filter: values})
    frame.to_csv(args.output, sep='\t', index=False)

if __name__ == '__main__':
    main()
