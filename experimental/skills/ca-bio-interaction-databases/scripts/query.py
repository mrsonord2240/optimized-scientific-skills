#!/usr/bin/env python3
"""Query interaction resources for a user-supplied gene list."""
import argparse
import os
from pathlib import Path
import pandas as pd
from interaction_clients import (aggregate_networks, biogrid_lt_physical,
                                 omnipath_interactions, signor_for_gene,
                                 string_network, string_resolve_ids)

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--resource', choices=['string', 'resolve', 'signor', 'biogrid', 'omnipath', 'aggregate'], required=True)
    p.add_argument('--genes', required=True, help='one gene symbol per line')
    p.add_argument('--out', required=True)
    p.add_argument('--caller')
    p.add_argument('--taxon', type=int, default=9606)
    p.add_argument('--threshold', type=int)
    p.add_argument('--network-type', choices=['functional', 'physical'])
    p.add_argument('--license', choices=['academic', 'commercial'], default='academic')
    a = p.parse_args()
    genes = [x.strip() for x in Path(a.genes).read_text(encoding='utf-8').splitlines() if x.strip()]
    if not genes:
        p.error('gene list is empty')
    if Path(a.out).exists():
        p.error('output exists; choose a fresh path')
    if a.resource not in ('string', 'resolve') and a.caller is not None:
        p.error('--caller applies to string and resolve only')
    if a.resource != 'string' and (a.threshold is not None or a.network_type is not None):
        p.error('--threshold and --network-type apply to string only')
    if a.resource != 'omnipath' and a.license != 'academic':
        p.error('--license commercial applies to omnipath only')
    if a.resource == 'aggregate' and (a.license != 'academic' or a.taxon != 9606):
        p.error('source aggregate supports its default academic human scope only')
    if a.resource in ('signor', 'omnipath') and a.taxon != 9606:
        p.error('source wrapper covers human SIGNOR/OmniPath usage')
    if a.resource == 'string':
        frame = string_network(genes, species=a.taxon, threshold=a.threshold if a.threshold is not None else 700,
                               network_type=a.network_type or 'functional', caller_identity=a.caller)
    elif a.resource == 'resolve':
        frame = string_resolve_ids(genes, species=a.taxon, caller_identity=a.caller)
    elif a.resource == 'signor':
        frame = pd.concat([signor_for_gene(g).assign(query_gene=g) for g in genes], ignore_index=True)
    elif a.resource == 'biogrid':
        key = os.environ.get('BIOGRID_API_KEY')
        if not key:
            p.error('BIOGRID_API_KEY environment variable is required')
        frame = pd.concat([biogrid_lt_physical(g, key, taxon=a.taxon).assign(query_gene=g) for g in genes], ignore_index=True)
    elif a.resource == 'omnipath':
        frame = omnipath_interactions(genes, license=a.license)
    else:
        graph = aggregate_networks(genes, biogrid_key=os.environ.get('BIOGRID_API_KEY'))
        frame = pd.DataFrame([{'gene_a': x, 'gene_b': y, 'sources': ','.join(sorted(v['sources'])),
                               'string_score': v['string_score']} for x, y, v in graph.edges(data=True)],
                             columns=['gene_a', 'gene_b', 'sources', 'string_score'])
    frame.to_csv(a.out, sep='\t', index=False)

if __name__ == '__main__':
    main()
