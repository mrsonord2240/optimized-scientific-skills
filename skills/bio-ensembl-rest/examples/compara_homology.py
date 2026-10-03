'''Ensembl Compara ortholog/paralog queries via REST; includes confidence handling.'''
# Reference: requests 2.31+, Ensembl REST release 116 | Verify API if version differs
import os
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'scripts'))
from ensembl_client import SLEEP, orthologs_by_symbol, paralogs_by_symbol


print('=== BRCA1 human -> all Compara orthologs ===')
all_orth = orthologs_by_symbol('human', 'BRCA1')
type_counts = {}
for h in all_orth:
    type_counts[h['type']] = type_counts.get(h['type'], 0) + 1
for t, c in sorted(type_counts.items(), key=lambda x: -x[1]):
    print(f'  {t:<28} {c}')
time.sleep(SLEEP)

print('\n=== BRCA1 human -> mouse, with confidence ===')
mouse_orth = orthologs_by_symbol('human', 'BRCA1', target='mouse')
for h in mouse_orth:
    print(f'  type={h["type"]}  confidence={h.get("confidence")}  '
          f'target={h["target"]["id"]}  '
          f'identity_target={h["target"].get("perc_id")}  '
          f'identity_query={h["source"].get("perc_id")}')
time.sleep(SLEEP)

print('\n=== BRCA1 within-species paralogs ===')
paralogs = paralogs_by_symbol('human', 'BRCA1')
if not paralogs:
    print('  none returned')
for p in paralogs[:5]:
    print(f'  paralog: {p["target"]["id"]}  type={p["type"]}  '
          f'taxonomic_level={p.get("taxonomy_level")}')
time.sleep(SLEEP)

print('\n=== Compara confidence semantics ===')
print('  type=ortholog_one2one: 1:1 ortholog (high-confidence single match)')
print('  type=ortholog_one2many: lineage-specific duplication in the target species')
print('  type=ortholog_many2many: ancestral duplication; multiple co-orthologs')
print('  type=within_species_paralog: in-paralog (post-speciation duplication)')
print('  confidence: often absent (None); when present it is 0 or 1. Treat missing as unknown.')
