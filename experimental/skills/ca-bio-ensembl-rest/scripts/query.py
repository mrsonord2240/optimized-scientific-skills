#!/usr/bin/env python3
"""Fetch one Ensembl REST JSON endpoint using the source retry client."""
import argparse
import json
from pathlib import Path
from ensembl_client import BASE, get_with_retry

def main():
    p = argparse.ArgumentParser()
    p.add_argument('path', help='endpoint path beginning with /')
    p.add_argument('output')
    p.add_argument('--base', default=BASE)
    p.add_argument('--param', action='append', default=[], help='query key=value; repeat as needed')
    args = p.parse_args()
    if not args.path.startswith('/') or args.path.startswith('//'):
        p.error('endpoint path must start with one slash')
    params = {}
    for item in args.param:
        if '=' not in item:
            p.error('parameters must be key=value')
        key, value = item.split('=', 1)
        params[key] = value
    response = get_with_retry(args.base.rstrip('/') + args.path, params=params)
    Path(args.output).write_text(json.dumps(response.json(), indent=2) + '\n', encoding='utf-8')

if __name__ == '__main__':
    main()
