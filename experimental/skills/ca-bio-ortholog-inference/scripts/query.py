#!/usr/bin/env python3
"""File-based entry points for the source ortholog clients."""
import argparse
import json
from pathlib import Path
from ortholog_clients import batch_compara, orthodb_search, orthodb_orthologs, oma_orthologs, ko_for_gene, genes_for_ko

parser = argparse.ArgumentParser()
sub = parser.add_subparsers(dest='task', required=True)
p = sub.add_parser('compara')
p.add_argument('--symbols', nargs='+', required=True)
p.add_argument('--source', default='human')
p.add_argument('--target', default='mouse')
p = sub.add_parser('orthodb-search')
p.add_argument('--query', required=True)
p.add_argument('--taxid', type=int, default=9606)
p = sub.add_parser('orthodb-members')
p.add_argument('--group', required=True)
p.add_argument('--taxid', type=int, default=10090)
p = sub.add_parser('oma')
p.add_argument('--protein', required=True)
p.add_argument('--rel-type', default=None)
p = sub.add_parser('kegg')
p.add_argument('--species', required=True)
p.add_argument('--gene', required=True)
for p in sub.choices.values():
    p.add_argument('--out', required=True)
a = parser.parse_args()
if a.task == 'compara':
    batch_compara(a.symbols, source=a.source, target=a.target).to_csv(a.out, index=False)
else:
    if a.task == 'orthodb-search':
        result = orthodb_search(a.query, a.taxid)
    elif a.task == 'orthodb-members':
        result = orthodb_orthologs(a.group, a.taxid)
    elif a.task == 'oma':
        result = oma_orthologs(a.protein, a.rel_type)
    else:
        ids = ko_for_gene(a.species, a.gene)
        result = {ko: genes_for_ko(ko) for ko in ids}
    Path(a.out).write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')
