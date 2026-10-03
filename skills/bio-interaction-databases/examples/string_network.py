'''STRING v12 network with per-channel inspection; functional-vs-physical filtering; multi-confidence comparison.'''
# Reference: requests 2.31+, pandas 2.2+, networkx 3.2+ | Verify API if version differs
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'scripts'))
import networkx as nx
import pandas as pd
from interaction_clients import string_network, string_resolve_ids


GENES = ['TP53', 'BRCA1', 'MDM2', 'ATM', 'CHEK2', 'CDK2', 'RB1', 'CDKN1A', 'BAX', 'BCL2']

print('=== Resolve gene symbols to STRING IDs ===')
print(string_resolve_ids(GENES)[['queryIndex', 'preferredName', 'stringId']].to_string(index=False))

print('\n=== Compare confidence tiers (number of edges retained) ===')
for thr in [400, 700, 900]:
    df = string_network(GENES, threshold=thr)
    print(f'  Threshold {thr}: {len(df)} edges')

print('\n=== Per-channel breakdown at threshold 700 ===')
df = string_network(GENES, threshold=700)
# Channels: experiments (escore), database (dscore), textmining (tscore),
# coexpression (ascore), neighborhood (nscore), fusion (fscore), cooccurrence (pscore)
df_show = df[['preferredName_A', 'preferredName_B', 'score',
              'escore', 'dscore', 'tscore', 'ascore']].head(8)
print(df_show.to_string(index=False))

print('\n=== Physical-only network (network_type=physical) ===')
# escore (experiments channel) mixes physical and functional experimental evidence, so it is
# not a physical filter; ask STRING for its physical network instead.
physical = string_network(GENES, threshold=700, network_type='physical')
print(f'  Functional network, combined score 700+: {len(df)} edges')
print(f'  Physical network, combined score 700+: {len(physical)} edges')
print('  Use the physical network for "physically interact" claims.')

print('\n=== Build NetworkX graph and rank by centrality ===')
g = nx.Graph()
for _, row in df.iterrows():
    g.add_edge(row['preferredName_A'], row['preferredName_B'],
               score=row['score'], escore=row['escore'])

print(f'  Nodes: {g.number_of_nodes()}; edges: {g.number_of_edges()}; '
      f'density: {nx.density(g):.3f}; components: {nx.number_connected_components(g)}')

deg = pd.DataFrame(sorted(g.degree(), key=lambda x: -x[1]), columns=['gene', 'degree'])
print('\n  Hubs (top 5 by degree):')
print(deg.head(5).to_string(index=False))
