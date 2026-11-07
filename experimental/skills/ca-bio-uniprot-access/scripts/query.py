#!/usr/bin/env python3
"""File-based entry points for the UniProt client."""
import argparse
import json
from pathlib import Path
import uniprot_client as client

parser = argparse.ArgumentParser()
sub = parser.add_subparsers(dest='task', required=True)
for name in ('entry', 'xrefs', 'isoforms'):
    sub.add_parser(name).add_argument('--accession', required=True)
for name in ('search', 'stream'):
    p = sub.add_parser(name)
    p.add_argument('--query', required=True)
    p.add_argument('--fields', required=True)
p = sub.add_parser('map')
p.add_argument('--ids', nargs='+', required=True)
p.add_argument('--from-db', default='Ensembl')
p.add_argument('--to-db', default='UniProtKB')
sub.add_parser('uniref').add_argument('--cluster', required=True)
sub.add_parser('proteome').add_argument('--id', required=True)
for p in sub.choices.values():
    p.add_argument('--out', required=True)
a = parser.parse_args()
if a.task in ('search', 'stream'):
    fn = client.search_tsv if a.task == 'search' else client.stream_tsv
    fn(a.query, a.fields.split(',')).to_csv(a.out, sep='\t', index=False)
elif a.task == 'proteome':
    client.download_proteome(a.id, a.out)
else:
    if a.task == 'entry':
        result = client.fetch_entry_json(a.accession)
    elif a.task == 'xrefs':
        result = client.xref_summary(a.accession)
    elif a.task == 'isoforms':
        result = client.list_isoforms(a.accession)
    elif a.task == 'map':
        result = client.map_ids(a.ids, from_db=a.from_db, to_db=a.to_db)
    else:
        result = client.uniref_cluster(a.cluster)
    Path(a.out).write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')
