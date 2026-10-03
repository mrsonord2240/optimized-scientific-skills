'''UniProt REST workflows: single entry parsing, bulk TSV search, /stream for >500, async ID mapping.'''
# Reference: requests 2.31+, pandas 2.2+, UniProt REST 2026_03 | Verify API if version differs
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'scripts'))
from uniprot_client import fetch_entry_json, map_ids, search_tsv, stream_tsv


print('=== Single entry parse (P04637 / TP53) ===')
for k, v in fetch_entry_json('P04637').items():
    print(f'  {k:<14} {v}')

print('\n=== Bulk search with fields= (human reviewed kinases) ===')
df, total = search_tsv(
    'organism_id:9606 AND reviewed:true AND keyword:KW-0418',
    fields=['accession', 'gene_primary', 'protein_name', 'length', 'xref_pdb'],
    size=500, with_total=True,
)
print(f'  {len(df)} rows returned of {total} matches (one page; use stream_tsv for all)')
print(df.head(5).to_string(index=False))

print('\n=== Stream for >500 results (all human Swiss-Prot reviewed) ===')
df_all = stream_tsv(
    'organism_id:9606 AND reviewed:true',
    fields=['accession', 'gene_primary', 'length'],
)
print(f'  All human Swiss-Prot: {len(df_all)}')

print('\n=== ID Mapping: Ensembl Gene -> UniProt ===')
mapping = map_ids(['ENSG00000141510', 'ENSG00000171862', 'ENSG00000139618'],
                   to_db='UniProtKB-Swiss-Prot')  # reviewed only; 'UniProtKB' adds TrEMBL rows
for r in mapping.get('results', []):
    print(f'  {r["from"]:<22} -> {r["to"]}')
for failed in mapping.get('failedIds', []):
    print(f'  {failed:<22} -> NOT MAPPED')
