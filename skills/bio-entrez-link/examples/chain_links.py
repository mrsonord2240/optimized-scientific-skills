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


print('=== Chain: Gene -> Protein (RefSeq) -> Structure ===')
tp53 = '7157'
proteins = link_one('gene', 'protein', tp53, linkname='gene_protein_refseq')
print(f'TP53 -> {len(proteins)} RefSeq proteins')
time.sleep(DELAY)

if proteins:
    protein_batch = ','.join(proteins[:10])
    h = Entrez.elink(dbfrom='protein', db='structure', id=protein_batch)
    r = Entrez.read(h); h.close()
    structures = []
    for ls in r:
        if ls['LinkSetDb']:
            structures.extend([l['Id'] for l in ls['LinkSetDb'][0]['Link']])
    print(f'  -> {len(structures)} structure UIDs in PDB-MMDB')
time.sleep(DELAY)

print('\n=== Large batch via history server (simulated) ===')
# 250 distinct gene UIDs -> two EPost chunks (200 + 50) that must be unioned
batch_genes = [str(i) for i in range(1, 251)]
print(f'Input: {len(batch_genes)} distinct gene UIDs')
we, qk = link_batch_via_history('gene', 'protein', batch_genes, linkname='gene_protein_refseq')
print(f'Got WebEnv (truncated): {we[:30] if we else "<none>"}...  QueryKey: {qk}')
h = Entrez.esearch(db='protein', term=f'#{qk}', WebEnv=we, usehistory='y', retmax=0)
n_linked = int(Entrez.read(h)['Count']); h.close()
# Linking from only the last 50-ID chunk gave ~200 proteins; the union covers all 250 genes.
print(f'History set #{qk} holds {n_linked} linked proteins')
assert n_linked > 1000, 'linked set is far too small -- earlier EPost chunks were dropped'
print('Downstream: Entrez.efetch(db=protein, WebEnv=we, query_key=qk, retstart=..., retmax=500) in batches')
