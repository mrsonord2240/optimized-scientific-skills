'''Chain ELink calls and demonstrate EPost + neighbor_history for batches >200 IDs.'''
# Reference: biopython 1.83+, entrez direct 21.0+ | Verify API if version differs
from Bio import Entrez
import time

Entrez.email = 'your.email@example.com'
DELAY = 0.34


def link_one(dbfrom, db, source_id, linkname=None):
    kwargs = {'dbfrom': dbfrom, 'db': db, 'id': source_id}
    if linkname:
        kwargs['linkname'] = linkname
    h = Entrez.elink(**kwargs)
    r = Entrez.read(h); h.close()
    if not r[0]['LinkSetDb']:
        return []
    return [link['Id'] for link in r[0]['LinkSetDb'][0]['Link']]


def link_batch_via_history(dbfrom, db, source_ids, linkname=None, chunk=200):
    '''EPost in chunks, then ELink with neighbor_history over ALL posted IDs.

    Each EPost creates its own QueryKey holding only that chunk (posting into an existing
    WebEnv does not merge), so linking from the last key alone silently drops every earlier
    chunk. The chunk keys are unioned with ESearch '#1 OR #2 ...' first. Chunk size 200 is
    conservative, not a limit: one EPost of 1,500 gene UIDs also succeeded (2026-09-21).
    '''
    webenv = None
    query_key = None
    keys = []
    for i in range(0, len(source_ids), chunk):
        kwargs = {'db': dbfrom, 'id': ','.join(source_ids[i:i+chunk])}
        if webenv:
            kwargs['WebEnv'] = webenv
        h = Entrez.epost(**kwargs)
        r = Entrez.read(h); h.close()
        webenv = r['WebEnv']
        query_key = r['QueryKey']
        keys.append(query_key)
        time.sleep(DELAY)

    if len(keys) > 1:
        h = Entrez.esearch(db=dbfrom, term=' OR '.join(f'#{k}' for k in keys),
                           WebEnv=webenv, usehistory='y', retmax=0)
        r = Entrez.read(h); h.close()
        query_key = r['QueryKey']
        time.sleep(DELAY)

    kwargs = {'dbfrom': dbfrom, 'db': db, 'cmd': 'neighbor_history',
              'WebEnv': webenv, 'query_key': query_key}
    if linkname:
        kwargs['linkname'] = linkname
    h = Entrez.elink(**kwargs)
    r = Entrez.read(h); h.close()
    if not r[0]['LinkSetDbHistory']:
        return None, None
    # WebEnv is top-level; QueryKey is per-LinkSetDbHistory entry.
    return r[0]['WebEnv'], r[0]['LinkSetDbHistory'][0]['QueryKey']


if __name__ == '__main__':
    import argparse
    import json
    import os
    from pathlib import Path
    p = argparse.ArgumentParser()
    p.add_argument('--email', required=True)
    p.add_argument('--source', required=True)
    p.add_argument('--target', required=True)
    p.add_argument('--ids-file', required=True)
    p.add_argument('--linkname')
    p.add_argument('--out', required=True)
    a = p.parse_args()
    ids = [x.strip() for x in Path(a.ids_file).read_text(encoding='utf-8').splitlines() if x.strip()]
    if not ids:
        p.error('input UID list is empty')
    if Path(a.out).exists():
        p.error('output exists; choose a fresh path')
    Entrez.email = a.email
    Entrez.api_key = os.environ.get('NCBI_API_KEY')
    webenv, key = link_batch_via_history(a.source, a.target, ids, a.linkname)
    Path(a.out).write_text(json.dumps({'source': a.source, 'target': a.target,
                                      'source_count': len(ids), 'webenv': webenv,
                                      'query_key': key}, indent=2) + '\n', encoding='utf-8')
