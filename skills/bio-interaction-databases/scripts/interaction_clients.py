"""Clients for interaction resources: STRING, BioGRID, SIGNOR, OmniPath, plus network helpers.

Usage:
    import sys; sys.path.insert(0, 'scripts')
    from interaction_clients import string_network, omnipath_interactions, aggregate_networks

Requires: requests 2.31+, pandas 2.2+, networkx 3.2+.
BioGRID needs a free API key (https://webservice.thebiogrid.org/); the others need none.
SIGNOR symbol lookup also calls the UniProt REST API (no key).

Every request has a timeout. Errors are re-raised without the request URL so a BioGRID
access key never appears in a message or traceback.
"""
import time
import warnings
from io import StringIO

import networkx as nx
import pandas as pd
import requests

STRING = 'https://version-12-5.string-db.org/api'   # pinned; use https://string-db.org/api to follow current
BIOGRID = 'https://webservice.thebiogrid.org/interactions/'
SIGNOR = 'https://signor.uniroma2.it/getData.php'
OMNI = 'https://omnipathdb.org'
UNIPROT = 'https://rest.uniprot.org/uniprotkb/search'
CALLER = 'bioskills-2026'   # default STRING caller_identity; pass your own app/user name
TIMEOUT = 60                # seconds per request
STRING_PAUSE = 1.0          # minimum seconds between STRING calls (one request at a time)

# BioGRID experimental systems counted as physical evidence.
PHYSICAL_LT_SYSTEMS = {
    'Affinity Capture-MS', 'Affinity Capture-Western', 'Affinity Capture-RNA',
    'Co-fractionation', 'Co-purification', 'Reconstituted Complex',
    'Co-crystal Structure', 'Two-hybrid', 'Far Western', 'FRET', 'PCA',
}

SIGNOR_COLUMNS = ['source', 'target', 'effect', 'mechanism', 'residue', 'pmid', 'score', 'signor_id']

_last_string_call = [0.0]
_omni_resources = {}


def _get(url, params, label, timeout=TIMEOUT):
    """GET with timeout; failures are re-raised without the URL (it may carry an access key)."""
    problem = None
    try:
        r = requests.get(url, params=params, timeout=timeout)
        r.raise_for_status()
        return r
    except requests.HTTPError as e:
        problem = f'{label}: HTTP {e.response.status_code} {e.response.reason}'
    except requests.RequestException as e:
        problem = f'{label}: {type(e).__name__}'
    raise requests.RequestException(problem)


def _string_get(path, params, timeout=TIMEOUT):
    wait = STRING_PAUSE - (time.monotonic() - _last_string_call[0])
    if wait > 0:
        time.sleep(wait)
    try:
        return _get(f'{STRING}/{path}', params, f'STRING {path}', timeout)
    finally:
        _last_string_call[0] = time.monotonic()


def string_network(genes, species=9606, threshold=700, network_type='functional',
                   caller_identity=None):
    """threshold: 150 (low), 400 (medium), 700 (high), 900 (highest).

    network_type: 'functional' (default) or 'physical' (STRING's physical-complex network;
    use it for "physically interact" claims, not escore).
    Columns: stringId_A/B, preferredName_A/B, ncbiTaxonId, score, nscore (neighborhood),
    fscore (fusion), pscore (cooccurrence), ascore (coexpression), escore (experiments),
    dscore (database), tscore (textmining). Scores are 0-1. Only edges among the query
    genes are returned.
    """
    params = {'identifiers': '%0d'.join(genes), 'species': species, 'required_score': threshold,
              'network_type': network_type, 'caller_identity': caller_identity or CALLER}
    r = _string_get('tsv/network', params)
    return pd.read_csv(StringIO(r.text), sep='\t')


def string_resolve_ids(genes, species=9606, caller_identity=None):
    """Columns: queryIndex (0-based position in `genes`), stringId, ncbiTaxonId, taxonName,
    preferredName, annotation."""
    params = {'identifiers': '%0d'.join(genes), 'species': species,
              'caller_identity': caller_identity or CALLER}
    r = _string_get('tsv/get_string_ids', params)
    return pd.read_csv(StringIO(r.text), sep='\t')


def biogrid_lt_physical(gene, api_key, taxon=9606):
    """Low-throughput physical interactions for one gene.

    One row per experiment (a pair with several papers or systems repeats), not per pair:
    use drop_duplicates on the sorted (gene_a, gene_b) pair for pair counts.
    A rejected or missing key raises RequestException with HTTP 401/403 (no URL in the text);
    an empty frame means no LT physical rows for that gene and taxon.
    """
    params = {'accesskey': api_key, 'format': 'json', 'searchNames': True, 'geneList': gene,
              'taxId': taxon, 'includeInteractors': True, 'max': 10000}
    r = _get(BIOGRID, params, 'BioGRID interactions')
    rows = []
    for v in r.json().values():
        if v['THROUGHPUT'] == 'Low Throughput' and v['EXPERIMENTAL_SYSTEM'] in PHYSICAL_LT_SYSTEMS:
            rows.append({'gene_a': v['OFFICIAL_SYMBOL_A'], 'gene_b': v['OFFICIAL_SYMBOL_B'],
                         'system': v['EXPERIMENTAL_SYSTEM'], 'pmid': v['PUBMED_ID']})
    return pd.DataFrame(rows, columns=['gene_a', 'gene_b', 'system', 'pmid'])


def uniprot_accession(gene_symbol, taxon=9606):
    """Reviewed UniProt accession for an exact gene symbol (e.g. TP53 -> P04637)."""
    params = {'query': f'gene_exact:{gene_symbol} AND organism_id:{taxon} AND reviewed:true',
              'fields': 'accession', 'format': 'tsv', 'size': 1}
    lines = _get(UNIPROT, params, 'UniProt search').text.strip().split('\n')
    if len(lines) < 2:
        raise ValueError(f'no reviewed UniProt entry for gene {gene_symbol!r} (taxon {taxon})')
    return lines[1].strip()


def signor_for_gene(gene_symbol, uniprot=None):
    """Signed, directed SIGNOR records that involve the gene's protein (either end).

    The symbol is resolved to a reviewed human UniProt accession (pass `uniprot` to skip).
    SIGNOR rows are per record (site/paper), so a pair can repeat. `effect` is SIGNOR's
    free text (e.g. 'up-regulates quantity by stabilization'); `score` is SIGNOR's 0-1 score.
    Returns an empty frame (with a warning) when SIGNOR answers 'No result found.'; raises
    ValueError when a non-empty answer has no parseable record rows (error page, format change).
    """
    acc = uniprot or uniprot_accession(gene_symbol)
    r = _get(SIGNOR, {'organism': 9606, 'id': acc}, 'SIGNOR getData')
    text = r.text.strip()
    if not text or text.startswith('No result found'):
        warnings.warn(f'SIGNOR returned no records for {gene_symbol} ({acc})')
        return pd.DataFrame(columns=SIGNOR_COLUMNS)
    rows = []
    for line in text.split('\n'):          # headerless, 29 tab-separated columns
        c = line.split('\t')
        if len(c) >= 28:
            rows.append({'source': c[0], 'target': c[4], 'effect': c[8], 'mechanism': c[9],
                         'residue': c[10], 'pmid': c[21], 'score': c[27], 'signor_id': c[26]})
    if not rows:
        raise ValueError(f'SIGNOR answer for {gene_symbol} ({acc}) is not parseable '
                         f'(no 28+ column rows): {text[:80]!r}')
    return pd.DataFrame(rows, columns=SIGNOR_COLUMNS)


def _omni_resource_purpose():
    """OmniPath resource name -> license purpose ('commercial', 'academic', 'non_profit', ...)."""
    if not _omni_resources:
        data = _get(f'{OMNI}/resources', {'format': 'json'}, 'OmniPath resources').json()
        for name, info in data.items():
            lic = info.get('license')
            _omni_resources[name] = lic.get('purpose') if isinstance(lic, dict) else None
    return _omni_resources


def _commercial_ok(token, purposes):
    name = token if token in purposes else token.split('_')[0]
    return purposes.get(name) == 'commercial'


def omnipath_interactions(genes, types=None, license='academic'):
    """types: post_translational, transcriptional, mirna_target, lncrna_target (None = default dataset).

    Returns each gene's whole neighbourhood (partners=), not only edges among `genes`.
    license='academic' returns everything. license='commercial' screens client-side
    (the server's license parameter does not filter): sources whose OmniPath license purpose
    is not 'commercial' (or unknown/composite) are removed from `sources`/`references` and rows
    left without any permitted source are dropped. This is a source-level screen only:
    is_directed/is_stimulation/consensus_* and curation_effort are computed server-side from all
    sources, and per-resource terms (NC, ND, share-alike) still need review.
    """
    if license not in ('academic', 'commercial'):
        raise ValueError("license must be 'academic' or 'commercial'")
    params = {'genesymbols': 1, 'fields': 'sources,references,curation_effort',
              'partners': ','.join(genes)}
    if types:
        params['types'] = types
    r = _get(f'{OMNI}/interactions', params, 'OmniPath interactions')
    df = pd.read_csv(StringIO(r.text), sep='\t')
    if license == 'commercial' and len(df):
        purposes = _omni_resource_purpose()
        keep, sources, refs = [], [], []
        for s, ref in zip(df['sources'], df['references'].fillna('')):
            ok = [t for t in str(s).split(';') if _commercial_ok(t, purposes)]
            keep.append(bool(ok))
            sources.append(';'.join(ok))
            refs.append(';'.join(x for x in str(ref).split(';') if x
                                 and _commercial_ok(x.split(':')[0], purposes)))
        df = df.assign(sources=sources, references=refs)[keep].reset_index(drop=True)
    return df


def _add_edge(g, a, b, source, string_score=None):
    if a == b:      # self-interaction: not a network edge
        return
    edge = g.get_edge_data(a, b, default={'sources': set(), 'string_score': None})
    edge['sources'].add(source)
    if string_score is not None:
        edge['string_score'] = max(edge['string_score'] or 0, string_score)
    g.add_edge(a, b, **edge)


def aggregate_networks(genes, biogrid_key=None):
    """Undirected union of edges AMONG the query genes (STRING >=700, OmniPath
    post-translational, optional BioGRID LT physical).

    Neighbour-expanding sources (OmniPath, BioGRID) are restricted to pairs with both ends in
    `genes`, so the sources cover the same scope. Each edge carries a `sources` set and
    `string_score` (STRING's 0-1 combined score, None when STRING did not report the edge;
    it is not a cross-resource confidence). Self-interactions (homodimers: a gene with itself, reported by
    BioGRID) are dropped so density stays in [0, 1]; get them from biogrid_lt_physical. For directed/signed resources keep a DiGraph
    instead (see examples/interaction_query.py).
    """
    wanted = set(genes)
    g = nx.Graph()
    for _, row in string_network(genes, threshold=700).iterrows():
        a, b = sorted([row['preferredName_A'], row['preferredName_B']])
        _add_edge(g, a, b, 'STRING', row['score'])
    for _, row in omnipath_interactions(genes, types='post_translational').iterrows():
        if row['source_genesymbol'] in wanted and row['target_genesymbol'] in wanted:
            a, b = sorted([row['source_genesymbol'], row['target_genesymbol']])
            _add_edge(g, a, b, 'OmniPath')
    if biogrid_key:
        for gene in genes:
            for _, row in biogrid_lt_physical(gene, biogrid_key).iterrows():
                if row['gene_a'] in wanted and row['gene_b'] in wanted:
                    a, b = sorted([row['gene_a'], row['gene_b']])
                    _add_edge(g, a, b, 'BioGRID-LT-physical')
    return g


def summary(g):
    return {
        'nodes': g.number_of_nodes(),
        'edges': g.number_of_edges(),
        'density': nx.density(g),
        'components': nx.number_connected_components(g),
        'mean_degree': sum(dict(g.degree()).values()) / max(g.number_of_nodes(), 1),
    }


def multi_source_edges(g, min_sources=2):
    """Edges reported by several resources (sources overlap; not a calibrated confidence)."""
    return [(a, b, d) for a, b, d in g.edges(data=True) if len(d['sources']) >= min_sources]
