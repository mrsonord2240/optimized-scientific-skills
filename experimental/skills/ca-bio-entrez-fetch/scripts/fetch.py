#!/usr/bin/env python3
"""Retrieve one EFetch or ESummary response to a fresh file."""
import argparse
import os
from pathlib import Path
from Bio import Entrez

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--email', required=True)
    p.add_argument('--db', required=True)
    p.add_argument('--id', required=True)
    p.add_argument('--out', required=True)
    p.add_argument('--rettype')
    p.add_argument('--retmode', default='text')
    p.add_argument('--summary', action='store_true')
    a = p.parse_args()
    if Path(a.out).exists():
        p.error('output exists; choose a fresh path')
    Entrez.email = a.email
    Entrez.api_key = os.environ.get('NCBI_API_KEY')
    if a.summary:
        handle = Entrez.esummary(db=a.db, id=a.id)
    else:
        kw = {'db': a.db, 'id': a.id, 'retmode': a.retmode}
        if a.rettype:
            kw['rettype'] = a.rettype
        handle = Entrez.efetch(**kw)
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
