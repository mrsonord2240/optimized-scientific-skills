"""Clients for pre-computed ortholog resources: Ensembl Compara, OrthoDB, OMA, KEGG Orthology.

Usage:
    import sys; sys.path.insert(0, 'scripts')
    from ortholog_clients import compara_orthologs, batch_compara, orthodb_groups, oma_orthologs

Requires: requests 2.31+, pandas 2.2+. No API keys. Ensembl allows 15 req/sec (55K/hour):
every helper goes through get_with_retry (timeout, retry on 429/5xx/connection errors) and
batch_compara sleeps SLEEP between calls.
"""
import time

import pandas as pd
import requests

ENSEMBL = 'https://rest.ensembl.org'
ORTHODB = 'https://data.orthodb.org/v12'
OMA = 'https://omabrowser.org/api'
KEGG = 'https://rest.kegg.jp'
HEADERS = {'Accept': 'application/json'}
SLEEP = 0.07   # Ensembl 15 req/sec ceiling
TIMEOUT = (10, 45)   # (connect, read) seconds; Ensembl homology calls can hang without one
_SESSION = requests.Session()


def get_with_retry(url, params=None, max_retries=4, headers=None, allow_404=False):
    """GET through one session with a timeout. Retries 429 (honoring Retry-After), 5xx and
    connection/timeout errors with backoff; raises on 4xx and when retries are exhausted.
    allow_404=True returns the 404 response instead of raising."""
    hdrs = headers if headers is not None else (HEADERS if url.startswith(ENSEMBL) else {})
    last = None
    for attempt in range(max_retries):
        try:
            r = _SESSION.get(url, params=params, headers=hdrs, timeout=TIMEOUT)
        except (requests.ConnectionError, requests.Timeout) as e:
            last = e
            time.sleep(2 ** attempt)
            continue
        if r.status_code == 429:
            last = requests.HTTPError(f'429 from {url}', response=r)
            time.sleep(int(r.headers.get('Retry-After', '5')))
            continue
        if r.status_code >= 500:
            last = requests.HTTPError(f'{r.status_code} from {url}', response=r)
            time.sleep(2 ** attempt)
            continue
        if r.status_code == 404 and allow_404:
            return r
        r.raise_for_status()
        return r
    raise requests.RequestException(f'{max_retries} attempts failed for {url}: {last}')


# ---- Ensembl Compara ------------------------------------------------------------------

def resolve_symbol(species, symbol):
    """Resolve a symbol to an Ensembl Gene ID first -- symbols are unstable (MARCH1 -> MARCHF1)."""
    return get_with_retry(f'{ENSEMBL}/lookup/symbol/{species}/{symbol}').json()['id']


def compara_orthologs(species, symbol, target_species=None):
    """Orthologs of a gene by symbol. `type` is ortholog_one2one / one2many / many2many.
    Live Compara records carry no `confidence` key (verified 2026-10-03), so `confidence`
    is None unless a future release adds it; use `type`, `taxonomy_level` and the
    identity percentages (`source_pid`, `target_pid`) as the quality signals."""
    params = {'type': 'orthologues'}
    if target_species:
        params['target_species'] = target_species
    r = get_with_retry(f'{ENSEMBL}/homology/symbol/{species}/{symbol}', params=params)
    data = r.json()['data']
    if not data:
        return []
    return [{
        'source_id': h['source']['id'],
        'source_pid': h['source'].get('perc_id'),
        'target_species': h['target']['species'],
        'target_id': h['target']['id'],
        'target_pid': h['target'].get('perc_id'),
        'type': h['type'],
        'taxonomy_level': h.get('taxonomy_level'),
        'confidence': h.get('confidence'),
    } for h in data[0]['homologies']]


def compara_one2one(symbol, source='mouse', target='human'):
    """The 1:1 ortholog dict (type ortholog_one2one), or None."""
    hits = [o for o in compara_orthologs(source, symbol, target) if o['type'] == 'ortholog_one2one']
    return hits[0] if hits else None


def batch_compara(symbols, source='human', target='mouse', sleep=SLEEP):
    """Wide table of Compara orthologs for many symbols. Every input symbol appears with a
    `status`: 'found', 'no ortholog returned' (empty result) or 'request failed' (HTTP, timeout or
    connection error; message in the `error` column), so a partial batch cannot read as complete."""
    rows = []
    for sym in symbols:
        try:
            found = compara_orthologs(source, sym, target)
        except requests.RequestException as e:
            rows.append({'symbol': sym, 'status': 'request failed', 'error': str(e)})
        else:
            if found:
                rows.extend({'symbol': sym, 'status': 'found', **o} for o in found)
            else:
                rows.append({'symbol': sym, 'status': 'no ortholog returned'})
        time.sleep(sleep)
    cols = ['symbol', 'status', 'error', 'source_id', 'source_pid', 'target_species', 'target_id',
            'target_pid', 'type', 'taxonomy_level', 'confidence']
    return pd.DataFrame(rows, columns=cols)


# ---- OrthoDB v12 ----------------------------------------------------------------------

def orthodb_search(query, species_taxid=9606):
    """OrthoDB /search is FULL TEXT over group descriptions, not a gene-symbol lookup:
    'TP53' ranks phosphoglycerate mutase (TIGAR) and TP53-target groups first. Returns dicts
    (id, name, level_name, gene_count) in relevance order; verify a group before use
    (see orthodb_group_has_gene)."""
    r = get_with_retry(f'{ORTHODB}/search', params={'query': query, 'species': species_taxid})
    body = r.json()
    return [{'id': g.get('id'), 'name': g.get('name'), 'level_name': g.get('level_name'),
             'gene_count': g.get('gene_count')} for g in (body.get('bigdata') or [])]


def orthodb_groups(symbol, species_taxid=9606):
    """Group IDs from the full-text search, relevance order (not guaranteed to be the
    gene's own family)."""
    r = get_with_retry(f'{ORTHODB}/search', params={'query': symbol, 'species': species_taxid})
    return r.json().get('data') or []


def orthodb_orthologs(og_id, species_taxid=10090):
    """Gene IDs (symbols) of group members in the target species/level (10090 = mouse).
    /orthologs returns groups -> genes -> gene_id.id; `data` is null (-> []) when the group
    has no member at that level. Raises on HTTP errors."""
    r = get_with_retry(f'{ORTHODB}/orthologs', params={'id': og_id, 'species': species_taxid})
    ids = []
    for grp in r.json().get('data') or []:
        for g in grp.get('genes') or []:
            gid = (g.get('gene_id') or {}).get('id')
            if gid:
                ids.append(gid)
    return ids


def orthodb_group_has_gene(og_id, symbol, species_taxid=9606):
    """True if `symbol` (case-insensitive) is a member of the group among `species_taxid` genes."""
    return symbol.upper() in {i.upper() for i in orthodb_orthologs(og_id, species_taxid)}


# ---- OMA ------------------------------------------------------------------------------

def oma_orthologs(uniprot_acc, rel_type=None):
    """Ortholog dicts for a UniProt or OMA ID; [] on 404 or when OMA has no orthologs.
    Items carry omaid, canonicalid, species.taxon_id, species.code. rel_type='1:1' (or
    '1:n', 'm:1', 'm:n') filters server-side and is the most reliable call. An empty list
    can mean 'no OMA orthologs' (BRCA1 P38398), not a failure; failures raise."""
    r = get_with_retry(f'{OMA}/protein/{uniprot_acc}/orthologs/',
                       params={'rel_type': rel_type} if rel_type else None, allow_404=True)
    return [] if r.status_code == 404 else r.json()


def oma_hog_for_protein(oma_or_uniprot_id):
    return get_with_retry(f'{OMA}/protein/{oma_or_uniprot_id}/').json().get('oma_hog_id')


def oma_hog_members(hog_id, level=None):
    return get_with_retry(f'{OMA}/hog/{hog_id}/', params={'level': level} if level else None).json()


# ---- KEGG Orthology (plain-text TSV, not JSON) ----------------------------------------

def ko_for_gene(species, gene):
    """species: KEGG code (hsa=human, mmu=mouse, dme=fly). gene: NCBI Gene ID or KEGG locus."""
    r = get_with_retry(f'{KEGG}/link/ko/{species}:{gene}')
    return [line.split('\t')[1].replace('ko:', '') for line in r.text.strip().split('\n') if line]


def genes_for_ko(ko_id):
    """All KEGG genes annotated with this KO."""
    r = get_with_retry(f'{KEGG}/link/genes/{ko_id}')
    return [line.split('\t')[1] for line in r.text.strip().split('\n') if line]


def ko_info(ko_id):
    """Description + pathway list for a KO."""
    r = get_with_retry(f'{KEGG}/get/{ko_id}')
    return r.text
