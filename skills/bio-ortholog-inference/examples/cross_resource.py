'''Cross-resource ortholog consensus: query Ensembl Compara, OMA REST, and OrthoDB; surface disagreement.'''
# Reference: requests 2.31+ | Verify API if version differs
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'scripts'))
import time

from ortholog_clients import (compara_orthologs, oma_orthologs, orthodb_group_has_gene,
                              orthodb_orthologs, orthodb_search)

SLEEP = 0.1
SYMBOL, UNIPROT = 'TP53', 'P04637'   # OMA has no orthologs for BRCA1 P38398 (200 with [])


def run(label, fn):
    """Isolate each resource: report unavailable (raised) separately from empty (no call)."""
    try:
        out = fn()
    except Exception as e:  # network, 5xx after retries, schema surprise
        print(f'  UNAVAILABLE ({type(e).__name__}): {e}')
        return None
    if not out:
        print('  empty: the resource returned no call')
    return out


print(f'=== Compara: {SYMBOL} human -> mouse ===')
hits = run('compara', lambda: [o['target_id'] for o in compara_orthologs('human', SYMBOL, 'mouse')
                               if o['type'].startswith('ortholog')])
if hits:
    print(f'  {len(hits)} hits: {hits}')
time.sleep(SLEEP)

print(f'\n=== OMA: {SYMBOL} (UniProt {UNIPROT}) human -> mouse, 1:1 ===')
oma = run('oma', lambda: oma_orthologs(UNIPROT, rel_type='1:1'))
if oma:
    mouse = [o.get('canonicalid') or o.get('omaid') for o in oma if o.get('species', {}).get('code') == 'MOUSE']
    print(f'  OMA 1:1 orthologs: {len(oma)}; mouse: {mouse}')
time.sleep(SLEEP)

print(f'\n=== OrthoDB: {SYMBOL} (full-text search, group verified by human member) ===')


def orthodb_leg():
    # /search is full text: TP53 returns TIGAR/TP53-target groups first. Search by protein
    # name, keep groups named exactly that, then keep one that contains the human gene.
    for g in orthodb_search('tumor protein p53', species_taxid=9606):
        if g['name'].lower() != 'tumor protein p53':
            continue
        time.sleep(SLEEP)
        if orthodb_group_has_gene(g['id'], SYMBOL, 9606):
            return g, sorted(set(orthodb_orthologs(g['id'], 10090)))
    return None


found = run('orthodb', orthodb_leg)
if found:
    g, mouse_members = found
    print(f'  Verified group {g["id"]} ({g["name"]}, {g["level_name"]}); mouse members: {mouse_members[:5]}')

print('\n=== Interpretation ===')
print('Disagreement across resources is informative -- not error.')
print('Compara is tree-reconciled (best for vertebrates).')
print('OMA is strict (high precision, lower recall).')
print('OrthoDB is broad (broad coverage, more ambiguous calls).')
print('For publication-grade calls, intersect 2+ resources; inspect disagreements case by case.')
