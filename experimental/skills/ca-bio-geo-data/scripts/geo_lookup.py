#!/usr/bin/env python3
"""Retrieve GEO search, PubMed-linked records, or a GSE's source-pattern SRA runs."""
import argparse
import json
import os
from pathlib import Path
import time

def main():
    p = argparse.ArgumentParser()
    p.add_argument('operation', choices=['search', 'pubmed', 'sra'])
    p.add_argument('--email')
    p.add_argument('--term')
    p.add_argument('--pmid')
    p.add_argument('--gse')
    p.add_argument('--retmax', type=int, default=20)
    p.add_argument('--out', required=True)
    a = p.parse_args()
    if Path(a.out).exists():
        p.error('output exists; choose a fresh path')
    if a.operation == 'sra':
        if not a.gse:
            p.error('--gse is required')
        from pysradb import SRAweb
        db = SRAweb()
        studies = db.gse_to_srp(a.gse)
        if studies.empty:
            Path(a.out).write_text('run_accession\n', encoding='utf-8')
            return
        runs = db.srp_to_srr(studies['study_accession'].iloc[0])
        runs.to_csv(a.out, sep='\t', index=False)
        print('Selected first SRP; inspect additional study accessions:', studies['study_accession'].tolist())
        return
    if not a.email:
        p.error('--email is required for Entrez')
    from Bio import Entrez
    Entrez.email = a.email
    Entrez.api_key = os.environ.get('NCBI_API_KEY')
    def read(handle):
        try:
            return Entrez.read(handle)
        finally:
            handle.close()
    if a.operation == 'search':
        if not a.term or a.retmax < 0:
            p.error('search needs --term and nonnegative --retmax')
        result = read(Entrez.esearch(db='gds', term=a.term, retmax=a.retmax))
        ids = result['IdList']
    else:
        if not a.pmid:
            p.error('pubmed needs --pmid')
        result = read(Entrez.elink(dbfrom='pubmed', db='gds', id=a.pmid))
        ids = [link['Id'] for linkset in result for group in linkset.get('LinkSetDb', []) for link in group.get('Link', [])]
    if ids:
        time.sleep(0.1 if Entrez.api_key else 0.34)
    summaries = read(Entrez.esummary(db='gds', id=','.join(ids))) if ids else []
    Path(a.out).write_text(json.dumps({'query': result, 'summaries': summaries}, default=str, indent=2) + '\n', encoding='utf-8')

if __name__ == '__main__':
    main()
