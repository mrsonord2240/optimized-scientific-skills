#!/usr/bin/env python3
"""Resolve source UID labels and save a small ELink response as XML."""
import argparse
import os
from pathlib import Path
from Bio import Entrez
from check_source_ids import describe_ids

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--email', required=True)
    p.add_argument('--source', required=True)
    p.add_argument('--target')
    p.add_argument('--id', action='append', required=True)
    p.add_argument('--linkname')
    p.add_argument('--cmd', choices=['neighbor', 'neighbor_score', 'acheck', 'ncheck', 'lcheck', 'llinks', 'prlinks'], default='neighbor')
    p.add_argument('--out', required=True)
    a = p.parse_args()
    if Path(a.out).exists():
        p.error('output exists; choose a fresh path')
    if a.cmd in ('neighbor', 'neighbor_score') and not a.target:
        p.error('target is required for neighbor requests')
    Entrez.email = a.email
    Entrez.api_key = os.environ.get('NCBI_API_KEY')
    for uid, label in describe_ids(a.source, a.id).items():
        print(f'{a.source}:{uid}\t{label}')
    kwargs = {'dbfrom': a.source, 'id': a.id, 'cmd': a.cmd}
    if a.target:
        kwargs['db'] = a.target
    if a.linkname:
        kwargs['linkname'] = a.linkname
    handle = Entrez.elink(**kwargs)
    try:
        raw = handle.read()
    finally:
        handle.close()
    data = raw if isinstance(raw, bytes) else raw.encode('utf-8')
    if not data.strip() or data.lstrip().lower().startswith((b'<html', b'<!doctype html')):
        raise RuntimeError('empty or HTML response; no file written')
    Path(a.out).write_bytes(data)

if __name__ == '__main__':
    main()
