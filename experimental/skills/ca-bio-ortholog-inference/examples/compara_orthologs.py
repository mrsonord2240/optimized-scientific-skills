'''Ensembl Compara REST: single gene + batch, with timeouts, retries and ortholog-type semantics.'''
# Reference: requests 2.31+, Ensembl REST release 116 | Verify API if version differs
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'scripts'))
import requests

from ortholog_clients import batch_compara, compara_orthologs, resolve_symbol


print('=== Single-gene: BRCA1 human -> mouse ===')
for o in compara_orthologs('human', 'BRCA1', target_species='mouse'):
    print(f'  {o["target_id"]}  type={o["type"]}  level={o["taxonomy_level"]}  '
          f'identity={o["target_pid"]}%')

print('\n=== Symbol resolution example (MARCH1 was renamed in 2020) ===')
try:
    old = resolve_symbol('human', 'MARCH1')
    print(f'  MARCH1 -> {old}')
except requests.HTTPError:
    print('  MARCH1 not found (renamed to MARCHF1 by HGNC 2020)')
new = resolve_symbol('human', 'MARCHF1')
print(f'  MARCHF1 -> {new}')

print('\n=== Batch: 5 genes, human -> zebrafish ===')
df = batch_compara(['TP53', 'BRCA1', 'MYC', 'ATM', 'MDM2'], source='human', target='zebrafish')
one2one = df[df['type'] == 'ortholog_one2one']
print(one2one[['symbol', 'target_id', 'taxonomy_level', 'target_pid']].to_string(index=False))
print('\nOutcome per input symbol:')
for sym in ['TP53', 'BRCA1', 'MYC', 'ATM', 'MDM2']:
    sub = df[df['symbol'] == sym]
    status = sub['status'].iloc[0]
    n11 = int((sub['type'] == 'ortholog_one2one').sum())
    detail = {'found': f'{len(sub)} calls, {n11} 1:1', 'request failed': sub['error'].iloc[0]}.get(status, '')
    print(f'  {sym}: {status}' + (f' ({detail})' if detail else ''))
failed = int((df['status'] == 'request failed').sum())
empty = int((df['status'] == 'no ortholog returned').sum())
print(f'\n1:1 calls: {len(one2one)}; request failed: {failed}; no ortholog returned: {empty}'
      + ('  -> PARTIAL BATCH, do not read as complete' if failed or empty else ''))
