'''Fetch a large ESearch result set to a FASTA file via the NCBI history server (webenv/query_key), no UIDs re-sent.

Inputs: --email (required by NCBI), --term (Entrez query), --out (FASTA path); optional --db (default
nucleotide), --chunk (records per EFetch, default 500). NCBI_API_KEY in the environment raises the rate limit.
Usage: python history_fetch.py --email you@inst.edu --term 'Homo sapiens[ORGN] AND srcdb_refseq[PROP] AND biomol_mrna[PROP]' --out out.fasta
Tested: Biopython 1.88, live NCBI E-utilities 2026-09-21.
'''
import argparse
import os
import time
from Bio import Entrez


def history_fetch(term, out_path, db='nucleotide', chunk=500):
    h = Entrez.esearch(db=db, term=term, usehistory='y', retmax=0)
    r = Entrez.read(h); h.close()
    total = int(r['Count'])

    with open(out_path, 'w') as out:
        for start in range(0, total, chunk):
            h = Entrez.efetch(db=db, rettype='fasta', retmode='text',
                              retstart=start, retmax=chunk,
                              webenv=r['WebEnv'], query_key=r['QueryKey'])
            out.write(h.read()); h.close()
            time.sleep(0.1 if Entrez.api_key else 0.34)
    return total


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--email', required=True)
    ap.add_argument('--term', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--db', default='nucleotide')
    ap.add_argument('--chunk', type=int, default=500)
    a = ap.parse_args()
    Entrez.email = a.email
    Entrez.api_key = os.environ.get('NCBI_API_KEY')
    print(f'{history_fetch(a.term, a.out, a.db, a.chunk)} records requested -> {a.out}')
