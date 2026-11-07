#!/usr/bin/env python3
"""Search, inspect or spell-check Entrez and save the parsed response."""
import argparse
import json
import os
from pathlib import Path
import time
from Bio import Entrez

DATABASES = 'pubmed,pmc,nucleotide,protein,gene,sra,gds,bioproject,biosample,clinvar'

def main():
    p = argparse.ArgumentParser()
    p.add_argument('operation', choices=['search', 'info', 'spell', 'counts'])
    p.add_argument('--email', required=True)
    p.add_argument('--db')
    p.add_argument('--dbs', default=DATABASES)
    p.add_argument('--term')
    p.add_argument('--retmax', type=int, default=100)
    p.add_argument('--history', action='store_true')
    p.add_argument('--out', required=True)
    a = p.parse_args()
    if Path(a.out).exists():
        p.error('output exists; choose a fresh path')
    if a.operation in ('search', 'spell') and not a.db:
        p.error('--db is required')
    if a.operation != 'info' and (not a.term or not a.term.strip()):
        p.error('--term is required')
    if a.retmax < 0:
        p.error('--retmax must be nonnegative')
    if a.term and (len(a.term) > 2000 or a.term.upper().count(' OR ') > 100):
        p.error('query is too long; post large identifier sets instead')
    Entrez.email = a.email
    Entrez.api_key = os.environ.get('NCBI_API_KEY')
    def read(handle):
        try:
            return Entrez.read(handle)
        finally:
            handle.close()
    if a.operation == 'info':
        result = read(Entrez.einfo(db=a.db) if a.db else Entrez.einfo())
    elif a.operation == 'spell':
        result = read(Entrez.espell(db=a.db, term=a.term))
    elif a.operation == 'search':
        result = read(Entrez.esearch(db=a.db, term=a.term, retmax=a.retmax,
                                     usehistory='y' if a.history else 'n'))
    else:
        result = {}
        databases = [x.strip() for x in a.dbs.split(',') if x.strip()]
        if not databases:
            p.error('--dbs must name at least one database')
        for db in databases:
            result[db] = read(Entrez.esearch(db=db, term=a.term, retmax=0))
            time.sleep(0.1 if Entrez.api_key else 0.34)
    Path(a.out).write_text(json.dumps(result, default=str, indent=2) + '\n', encoding='utf-8')

if __name__ == '__main__':
    main()
