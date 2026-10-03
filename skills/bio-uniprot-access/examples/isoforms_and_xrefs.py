'''Iterate isoforms and cross-references for a UniProt entry; download reference proteome.'''
# Reference: requests 2.31+, UniProt REST 2026_03 | Verify API if version differs
import os
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'scripts'))
from uniprot_client import download_proteome, fetch_isoform_fasta, list_isoforms, xref_summary

DELAY = 0.05  # 200 req/sec tolerated


print('=== TP53 isoforms (P04637) ===')
isos = list_isoforms('P04637')
print(f'Found {len(isos)} isoforms:')
for iso in isos[:5]:
    canonical = '[canonical]' if iso['canonical'] else ''
    print(f'  {", ".join(iso["ids"]):<20} {iso["name"] or "":<30} {canonical}')

time.sleep(DELAY)

# Fetch a non-canonical isoform
if len(isos) > 1:
    second = isos[1]['ids'][0] if isos[1]['ids'] else None
    if second:
        fasta = fetch_isoform_fasta(second)
        print(f'\n=== Fetched isoform {second} FASTA (first 200 chars) ===')
        print(fasta[:200])
        time.sleep(DELAY)

print('\n=== Cross-references summary (P04637) ===')
xrefs = xref_summary('P04637')
for db in sorted(xrefs, key=lambda k: -len(xrefs[k]))[:10]:
    print(f'  {db:<20} {len(xrefs[db])} entries')

print('\n=== Download human reference proteome (commented out - large file) ===')
print('  download_proteome("UP000005640", "human_reference_proteome.fasta.gz")')
print('  # ~38 MB compressed; ~87 MB unpacked; ~147.5K entries (reviewed + TrEMBL)')
