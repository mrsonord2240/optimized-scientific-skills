"""Small Ensembl REST client: lookup, sequence, overlap, VEP, Compara, LD.

Usage:
    import sys; sys.path.insert(0, 'scripts')
    from ensembl_client import symbol_to_id, vep_region, orthologs_by_symbol

Every call goes through get_with_retry (30 s timeout; retries 429 with Retry-After,
5xx and timeouts; raises EnsemblError after 3 attempts or on a non-JSON body).
Pass base=ARCHIVE_E110 (or any eNNN archive host) for release pinning;
archive hosts answer with a redirect that requests follows automatically.
Requires: requests 2.31+. No API key.
"""
import time

import requests

BASE = 'https://rest.ensembl.org'
ARCHIVE_E110 = 'https://e110.rest.ensembl.org'   # redirects to the monthly archive host
GRCH37 = 'https://grch37.rest.ensembl.org'
HEADERS = {'Accept': 'application/json'}
SLEEP = 0.07   # 15 req/sec ceiling


class EnsemblError(RuntimeError):
    """Retries exhausted, or the server answered with a non-JSON body."""


def get_with_retry(url, params=None, max_retries=3, headers=None, timeout=30):
    """GET JSON. Retries 429 (honors Retry-After), 5xx and timeouts with backoff.

    4xx raises requests.HTTPError carrying the server's error message.
    Raises EnsemblError when retries run out or the body is not JSON (a retired
    archive can answer 200 with an HTML page).
    """
    last = 'no attempt made'
    for attempt in range(max_retries):
        try:
            r = requests.get(url, params=params, headers=headers or HEADERS, timeout=timeout)
        except (requests.Timeout, requests.ConnectionError) as e:
            last = type(e).__name__
            time.sleep(2 ** attempt)
            continue
        if r.status_code == 429:
            last = '429'
            try:
                wait = float(r.headers.get('Retry-After', '5'))
            except ValueError:
                wait = 5.0
            time.sleep(wait)
            continue
        if r.status_code >= 500:
            last = str(r.status_code)
            time.sleep(2 ** attempt)
            continue
        if r.status_code >= 400:
            raise requests.HTTPError(f'{r.status_code} {r.reason} for {r.url}: {r.text[:200]}',
                                     response=r)
        if 'json' not in r.headers.get('Content-Type', ''):
            raise EnsemblError(f'non-JSON response ({r.headers.get("Content-Type")}) from {r.url}')
        return r
    raise EnsemblError(f'{max_retries} attempts exhausted (last: {last}) for {url}')


def symbol_to_id(species, symbol, base=BASE):
    """Resolve a symbol to the canonical Ensembl gene record (use ['id'] downstream)."""
    return get_with_retry(f'{base}/lookup/symbol/{species}/{symbol}').json()


def gene_info(ensembl_id, base=BASE):
    """Gene record with transcripts and exons (expand=1)."""
    return get_with_retry(f'{base}/lookup/id/{ensembl_id}', params={'expand': 1}).json()


def sequence_for_id(ensembl_id, seq_type='protein', base=BASE, multiple_sequences=False):
    """seq_type: cdna, cds, protein, genomic.

    Use a transcript/translation ID (ENST/ENSP) for cdna, cds or protein. A gene ID
    (ENSG) with a non-genomic type needs multiple_sequences=True and then returns a
    list with one record per transcript instead of a single dict.
    """
    params = {'type': seq_type}
    if multiple_sequences:
        params['multiple_sequences'] = 1
    return get_with_retry(f'{base}/sequence/id/{ensembl_id}', params=params).json()


def genes_in_region(species, region, base=BASE, feature='gene'):
    """region as 'chr:start-end', e.g. '17:43000000-43200000'."""
    return get_with_retry(f'{base}/overlap/region/{species}/{region}',
                          params={'feature': feature}).json()


def vep_region(species, region, allele, base=BASE):
    """region format chr:start-end:strand; allele is the alt base(s)."""
    return get_with_retry(f'{base}/vep/{species}/region/{region}/{allele}').json()


def vep_hgvs(species, hgvs, base=BASE):
    return get_with_retry(f'{base}/vep/{species}/hgvs/{hgvs}').json()


def vep_id(species, variant_id, base=BASE):
    return get_with_retry(f'{base}/vep/{species}/id/{variant_id}').json()


def _homologies(species, symbol, params, base):
    r = get_with_retry(f'{base}/homology/symbol/{species}/{symbol}', params=params)
    data = r.json().get('data', [])
    return data[0].get('homologies', []) if data else []


def orthologs_by_symbol(species, symbol, target=None, base=BASE):
    params = {'type': 'orthologues'}
    if target:
        params['target_species'] = target
    return _homologies(species, symbol, params, base)


def paralogs_by_symbol(species, symbol, base=BASE):
    return _homologies(species, symbol, {'type': 'paralogues'}, base)


def ld_pairwise(species, var1, var2, population='1000GENOMES:phase_3:CEU', base=BASE):
    return get_with_retry(f'{base}/ld/{species}/pairwise/{var1}/{var2}',
                          params={'population_name': population}).json()


def batch_symbols(species, symbols, base=BASE):
    """Resolve many symbols; errors are recorded per symbol instead of raised."""
    out = {}
    for sym in symbols:
        try:
            out[sym] = symbol_to_id(species, sym, base=base)
        except (requests.RequestException, EnsemblError) as e:
            out[sym] = {'error': str(e)}
        time.sleep(SLEEP)
    return out
