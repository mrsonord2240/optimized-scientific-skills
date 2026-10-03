'''Ensembl REST lookup, sequence, overlap; demonstrates symbol -> ID resolution and archive pinning.'''
# Reference: requests 2.31+, Ensembl REST release 116 | Verify API if version differs
import os
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'scripts'))
from ensembl_client import (ARCHIVE_E110 as ARCHIVE, SLEEP, gene_info, genes_in_region,
                            sequence_for_id, symbol_to_id)


print('=== Symbol resolution (live release) ===')
info = symbol_to_id('human', 'BRCA1')
print(f'  Ensembl Gene ID: {info["id"]}')
print(f'  Biotype:         {info["biotype"]}')
print(f'  Location:        chr{info["seq_region_name"]}:{info["start"]}-{info["end"]} (strand {info["strand"]})')
time.sleep(SLEEP)

print('\n=== Same query against archive (release 110) -- for reproducibility ===')
info_pinned = symbol_to_id('human', 'BRCA1', base=ARCHIVE)
print(f'  Ensembl Gene ID (e110): {info_pinned["id"]}')
if info_pinned["id"] != info["id"]:
    print(f'  WARNING: live and archive Gene IDs differ -- gene model has changed')
time.sleep(SLEEP)

print('\n=== Full gene info (transcripts + exons via ?expand=1) ===')
detail = gene_info(info['id'])
print(f'  Transcripts: {len(detail.get("Transcript", []))}')
for tx in detail.get('Transcript', [])[:3]:
    print(f'    {tx["id"]:<22} {tx["biotype"]:<22} {len(tx.get("Exon", []))} exons')
time.sleep(SLEEP)

print('\n=== Protein sequence ===')
prot = sequence_for_id('ENSP00000269305', seq_type='protein')  # translation ID; a gene ID needs multiple_sequences=True
print(f'  TP53 protein {prot["id"]}: {len(prot["seq"])} aa  (first 60: {prot["seq"][:60]}...)')
time.sleep(SLEEP)

print('\n=== Genes in interval chr17:43000000-43200000 ===')
for g in genes_in_region('human', '17:43000000-43200000'):
    sym = g.get('external_name', '?')
    print(f'  {sym:<12} {g["id"]} {g["biotype"]:<22} {g["start"]}-{g["end"]}')
