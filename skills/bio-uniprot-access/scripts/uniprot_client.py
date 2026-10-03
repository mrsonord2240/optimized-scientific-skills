"""UniProt REST client (rest.uniprot.org): entries, search/stream TSV, ID mapping, isoforms, proteomes.

Usage:
    import sys; sys.path.insert(0, 'scripts')
    from uniprot_client import fetch_entry_json, search_tsv, stream_tsv, map_ids

Requires: requests 2.31+, pandas 2.2+. No API key.
Every request goes through one helper with a timeout, bounded retries on 429/5xx, and Retry-After handling.
"""
import time
import warnings
from io import StringIO

import pandas as pd
import requests

BASE = 'https://rest.uniprot.org'
TIMEOUT = 60  # seconds per request
RETRIES = 4


def _request(method, url, **kw):
    """requests wrapper: timeout, retry on 429/5xx (honours Retry-After), raise on other errors."""
    kw.setdefault('timeout', TIMEOUT)
    for attempt in range(RETRIES + 1):
        try:
            r = requests.request(method, url, **kw)
        except (requests.ConnectionError, requests.Timeout):
            if attempt == RETRIES:
                raise
            time.sleep(2 ** attempt)
            continue
        if r.status_code in (429, 500, 502, 503, 504) and attempt < RETRIES:
            try:
                wait = float(r.headers.get('Retry-After', ''))
            except ValueError:
                wait = 2 ** attempt
            r.close()
            time.sleep(min(wait, 60))
            continue
        r.raise_for_status()
        return r


def fetch_entry_json(accession):
    """One entry as a flat dict, using defensive .get() chains over the nested JSON schema.

    Secondary accessions redirect to the primary entry. Deleted (inactive) accessions return HTTP 200
    with entryType 'Inactive' and no sequence; those raise ValueError naming the reason.
    """
    e = _request('GET', f'{BASE}/uniprotkb/{accession}.json').json()
    if e.get('entryType') == 'Inactive':
        reason = e.get('inactiveReason', {})
        raise ValueError(f"{accession} is an inactive UniProtKB entry: {reason.get('inactiveReasonType')} "
                         f"({reason.get('deletedReason') or reason.get('mergeDemergeTo')})")
    xrefs = e.get('uniProtKBCrossReferences', [])
    pdb = [x['id'] for x in xrefs if x['database'] == 'PDB']
    return {
        'accession': e['primaryAccession'],
        'entry_name': e.get('uniProtkbId'),
        'reviewed': e.get('entryType', '').startswith('UniProtKB reviewed'),
        'protein_name': e.get('proteinDescription', {}).get('recommendedName', {}).get('fullName', {}).get('value'),
        'gene_primary': (e.get('genes') or [{}])[0].get('geneName', {}).get('value'),
        'sequence': e['sequence']['value'],
        'length': e['sequence']['length'],
        'pdb_count': len(pdb),
        'pdb_ids': pdb,
        'alphafold_id': next((x['id'] for x in xrefs if x['database'] == 'AlphaFoldDB'), None),
    }


def search_tsv(query, fields, size=500, with_total=False):
    """/search with explicit fields; first page only (max 500). Use stream_tsv for more.

    with_total=True returns (DataFrame, total) where total is the server's X-Total-Results.
    A UserWarning is emitted whenever the page is smaller than the total.
    """
    params = {'query': query, 'fields': ','.join(fields), 'format': 'tsv', 'size': size}
    r = _request('GET', f'{BASE}/uniprotkb/search', params=params)
    df = pd.read_csv(StringIO(r.text), sep='\t')
    total = int(r.headers.get('X-Total-Results', len(df)))
    if total > len(df):
        warnings.warn(f'search_tsv returned {len(df)} of {total} matches; use stream_tsv for the full set')
    return (df, total) if with_total else df


def stream_tsv(query, fields):
    """/stream has no 500-result cap; right endpoint for bulk pulls (response is held in memory)."""
    params = {'query': query, 'fields': ','.join(fields), 'format': 'tsv'}
    r = _request('GET', f'{BASE}/uniprotkb/stream', params=params)
    return pd.read_csv(StringIO(r.text), sep='\t')


def map_ids(ids, from_db='Ensembl', to_db='UniProtKB', timeout=600, poll_interval=3):
    """Async ID mapping: submit, poll with a hard timeout, fetch ALL results via /idmapping/stream.

    Returns {'results': [{'from', 'to'}, ...], 'failedIds': [...]}. failedIds holds the server's failures
    plus every input that produced no row. to_db='UniProtKB' includes TrEMBL rows (many per gene);
    use 'UniProtKB-Swiss-Prot' for reviewed entries only.
    """
    ids = list(ids)
    submit = _request('POST', f'{BASE}/idmapping/run',
                      data={'ids': ','.join(ids), 'from': from_db, 'to': to_db})
    job_id = submit.json()['jobId']

    start = time.monotonic()
    while True:
        status = _request('GET', f'{BASE}/idmapping/status/{job_id}').json()
        state = status.get('jobStatus')
        if state == 'FAILED':
            raise RuntimeError(f'ID mapping job {job_id} failed: {status}')
        if state not in ('RUNNING', 'NEW'):
            break  # FINISHED, or the results payload came back directly
        if time.monotonic() - start >= timeout:
            raise TimeoutError(f'ID mapping job {job_id} did not complete in {timeout}s')
        time.sleep(poll_interval)

    data = _request('GET', f'{BASE}/idmapping/stream/{job_id}').json()
    results = data.get('results', [])
    mapped = {r['from'] for r in results}
    failed = list(data.get('failedIds', []))
    failed += [i for i in ids if i not in mapped and i not in failed]
    return {'results': results, 'failedIds': failed}


def resolve_obsolete(accessions):
    """Update obsolete/merged accessions to current primary IDs (rows 'from' -> 'to' accession)."""
    return map_ids(accessions, from_db='UniProtKB_AC-ID', to_db='UniProtKB')


def list_isoforms(accession):
    """Read comments[type=ALTERNATIVE PRODUCTS] for isoform IDs and names ('ids' is a list)."""
    r = _request('GET', f'{BASE}/uniprotkb/{accession}.json')
    for comment in r.json().get('comments', []):
        if comment.get('commentType') == 'ALTERNATIVE PRODUCTS':
            return [{
                'name': iso.get('name', {}).get('value'),
                'ids': iso.get('isoformIds', []),
                'canonical': iso.get('isoformSequenceStatus') == 'Displayed',
            } for iso in comment.get('isoforms', [])]
    return []


def fetch_isoform_fasta(accession_with_suffix):
    """e.g. 'P04637-2'. The bare accession returns the canonical sequence only."""
    return _request('GET', f'{BASE}/uniprotkb/{accession_with_suffix}.fasta').text


def xref_summary(accession):
    """Cross-reference IDs grouped by database name."""
    r = _request('GET', f'{BASE}/uniprotkb/{accession}.json')
    by_db = {}
    for x in r.json().get('uniProtKBCrossReferences', []):
        by_db.setdefault(x['database'], []).append(x['id'])
    return by_db


def download_proteome(upid, out_path):
    """upid: UniProt Proteome ID, e.g. UP000005640 (human reference). Writes gzipped FASTA.

    There is no /proteomes/{upid}.fasta.gz route; the proteome is pulled through /uniprotkb/stream.
    """
    params = {'query': f'proteome:{upid}', 'format': 'fasta', 'compressed': 'true'}
    with _request('GET', f'{BASE}/uniprotkb/stream', params=params, stream=True) as r:
        with open(out_path, 'wb') as f:
            for chunk in r.iter_content(65536):
                f.write(chunk)
    return out_path


def uniref_cluster(uniref_id):
    """e.g. UniRef50_P04637 -- the UniRef50 cluster centered on P04637.

    identity is the numeric tier parsed from the id (50); representative is the first accession of the
    representative member, with its member id (entry name or UniParc id) alongside.
    """
    j = _request('GET', f'{BASE}/uniref/{uniref_id}.json').json()
    rep = j['representativeMember']
    accs = rep.get('accessions') or []
    return {
        'id': j['id'],
        'representative': accs[0] if accs else rep.get('memberId'),
        'representative_member_id': rep.get('memberId'),
        'member_count': j['memberCount'],
        'identity': int(j['id'].split('_')[0].replace('UniRef', '')),
    }
