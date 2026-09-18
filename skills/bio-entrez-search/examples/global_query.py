'''Cross-database counts via an ESearch loop (authoritative path).

Entrez.egquery() does not exist on Biopython >=1.85 (AttributeError -- the function was dropped from
the public API; Bio/Entrez/__init__.py's module docstring still lists it, but there is no `def
egquery`). Calling the underlying NCBI endpoint directly is not a working substitute either: as of
2026-09-17, https://eutils.ncbi.nlm.nih.gov/entrez/eutils/egquery.fcgi 301-redirects to
ext-http-eutils.linkerd.ncbi.nlm.nih.gov, an internal-only NCBI hostname that does not resolve
outside their network. This script demonstrates the failure explicitly, then uses the documented
fallback -- looping ESearch with retmax=0 over a curated database list -- as the primary path.
'''
# Reference: biopython 1.88 confirmed; egquery confirmed broken (Biopython + HTTP) 2026-09-17
from Bio import Entrez
import time

Entrez.email = 'your.email@example.com'

DELAY = 0.34  # 3 req/sec ceiling without API key
CURATED_DBS = ['pubmed', 'pmc', 'nucleotide', 'protein', 'gene', 'sra', 'gds', 'bioproject', 'biosample', 'clinvar']


def confirm_egquery_broken(term):
    '''Attempt Entrez.egquery() as SKILL.md's decision table used to document, and surface the
    failure instead of silently swallowing it.'''
    try:
        handle = Entrez.egquery(term=term)
        record = Entrez.read(handle); handle.close()
        return {r['DbName']: int(r['Count']) for r in record['eGQueryResult']}
    except AttributeError as e:
        print(f'Entrez.egquery() unavailable as documented: {e}')
        return None


def cross_db_counts(term, dbs=CURATED_DBS):
    '''Authoritative per-database counts via ESearch. Slower than a single EGQuery call would have
    been, but does not depend on the (currently broken) egquery.fcgi endpoint and never lags the
    per-database index.'''
    counts = {}
    for db in dbs:
        h = Entrez.esearch(db=db, term=term, retmax=0)
        r = Entrez.read(h); h.close()
        counts[db] = int(r['Count'])
        time.sleep(DELAY)
    return counts


if __name__ == '__main__':
    term = 'CRISPR'

    print('=== Step 1: confirm Entrez.egquery() is unavailable ===')
    confirm_egquery_broken(term)
    time.sleep(DELAY)

    print(f'\n=== Step 2: cross-database counts for {term!r} via ESearch loop ===')
    counts = cross_db_counts(term)
    for db in sorted(counts, key=counts.get, reverse=True):
        if counts[db] > 0:
            print(f'  {db:<15} {counts[db]:>10,}')
