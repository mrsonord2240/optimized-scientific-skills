'''Aggregate STRING + OmniPath + SIGNOR into a unified network with per-edge provenance and direction handling.'''
# Reference: requests 2.31+, pandas 2.2+, networkx 3.2+ | Verify API if version differs
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'scripts'))
import networkx as nx
import pandas as pd
from interaction_clients import omnipath_interactions, signor_for_gene, string_network


GENES = ['TP53', 'MDM2', 'BRCA1', 'ATM', 'CHEK2', 'CDK2', 'CDKN1A', 'RB1']

# Use DiGraph because SIGNOR and OmniPath are directional;
# STRING edges are undirected and added in both directions.
g = nx.DiGraph()


print('=== STRING (high confidence, combined score 700+) ===')
string_df = string_network(GENES, threshold=700)
for _, row in string_df.iterrows():
    a, b = row['preferredName_A'], row['preferredName_B']
    for src, tgt in [(a, b), (b, a)]:  # undirected -> two directed edges
        if g.has_edge(src, tgt):
            g[src][tgt]['sources'].add('STRING')
        else:
            g.add_edge(src, tgt, sources={'STRING'}, string_score=row['score'],
                       directional=False, signed_effect=None, mechanism=None)
print(f'  {len(string_df)} STRING edges added')


print('\n=== OmniPath (directional aggregate) ===')
omni_df = omnipath_interactions(GENES)
omni_local = omni_df[omni_df['source_genesymbol'].isin(GENES) & omni_df['target_genesymbol'].isin(GENES)]
for _, row in omni_local.iterrows():
    a, b = row['source_genesymbol'], row['target_genesymbol']
    if g.has_edge(a, b):
        g[a][b]['sources'].add('OmniPath')
        g[a][b]['directional'] = True
    else:
        g.add_edge(a, b, sources={'OmniPath'}, string_score=None,
                   directional=True, signed_effect=None, mechanism=None)
print(f'  {len(omni_local)} OmniPath directional edges added')


print('\n=== SIGNOR (signed, mechanism-typed) ===')
sig_count = 0
for gene in GENES:
    sig_df = signor_for_gene(gene)
    for _, row in sig_df.iterrows():
        a, b = row['source'], row['target']
        if a in GENES and b in GENES:
            if g.has_edge(a, b):
                g[a][b]['sources'].add('SIGNOR')
                g[a][b]['signed_effect'] = row['effect']
                g[a][b]['mechanism'] = row['mechanism']
            else:
                g.add_edge(a, b, sources={'SIGNOR'}, string_score=None,
                           directional=True, signed_effect=row['effect'],
                           mechanism=row['mechanism'])
            sig_count += 1
print(f'  {sig_count} SIGNOR records among the query genes added (records repeat per site/paper)')
if sig_count == 0:
    print('  WARNING: no SIGNOR records among these genes')


print('\n=== Aggregated network summary ===')
print(f'  Nodes: {g.number_of_nodes()}')
print(f'  Directed edges: {g.number_of_edges()}')
multi = [(a, b, d) for a, b, d in g.edges(data=True) if len(d['sources']) >= 2]
print(f'  Edges reported by 2+ resources: {len(multi)}')
signed = [(a, b, d) for a, b, d in g.edges(data=True) if d.get('signed_effect')]
print(f'  Signed (SIGNOR) edges: {len(signed)}')


print('\n=== Signed signaling edges with mechanism ===')
for a, b, d in signed[:8]:
    print(f'  {a} --{d["signed_effect"]}--> {b}  ({d["mechanism"]})')


print('\n=== Export edge list with provenance ===')
edge_rows = []
for a, b, d in g.edges(data=True):
    edge_rows.append({
        'source': a, 'target': b,
        'sources': ','.join(sorted(d['sources'])),
        'n_sources': len(d['sources']),
        'signed_effect': d.get('signed_effect'),
        'mechanism': d.get('mechanism'),
        'string_score': d['string_score'],   # STRING 0-1 combined score; empty for other resources
        'directional': d['directional'],
    })
pd.DataFrame(edge_rows).to_csv('aggregated_interactions.csv', index=False)
print(f'  Wrote aggregated_interactions.csv ({len(edge_rows)} edges)')
