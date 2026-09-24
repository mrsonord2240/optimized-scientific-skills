'''Compare retrieval strategies: direct EFetch, history server, and NCBI Datasets CLI, with cost/payload metrics.'''
# Verified: biopython 1.88, ncbi datasets cli 18.37.0 (2026-09-22) | Verify API if version differs
from Bio import Entrez, SeqIO
import os
import subprocess
import sys
import time

Entrez.email = 'your.email@example.com'
DELAY = 0.34
# Resolved from PATH by default; set $DATASETS if the CLI lives elsewhere.
DATASETS_BIN = os.environ.get('DATASETS', 'datasets')


def history_server_download(db, term, out_path, batch_size=500):
    h = Entrez.esearch(db=db, term=term, usehistory='y', retmax=0)
    s = Entrez.read(h); h.close()
    webenv, query_key, total = s['WebEnv'], s['QueryKey'], int(s['Count'])
    print(f'  history-server: {total} records')

    if total == 0:
        return 0

    t0 = time.time()
    # newline='' disables newline translation: text mode would otherwise turn
    # NCBI's LF endings into CRLF on Windows, changing the payload's bytes.
    with open(out_path, 'w', newline='') as out:
        for start in range(0, total, batch_size):
            h = Entrez.efetch(db=db, rettype='fasta', retmode='text',
                              retstart=start, retmax=batch_size,
                              webenv=webenv, query_key=query_key)
            out.write(h.read()); h.close()
            time.sleep(DELAY)
    elapsed = time.time() - t0
    return elapsed


def datasets_download(target, out_zip, taxon=False,
                      include='genome,protein,gff3', extra_args=()):
    '''The supported bulk path for genome data -- replaces assembly_summary.txt scraping.

    taxon=False resolves `target` as an accession; taxon=True treats it as a
    taxonomic name (add '--reference' via extra_args to restrict to reference genomes).
    '''
    cmd = [DATASETS_BIN, 'download', 'genome', 'taxon' if taxon else 'accession', target,
           '--include', include, '--filename', out_zip, *extra_args]
    t0 = time.time()
    try:
        subprocess.run(cmd, check=True)
    except FileNotFoundError:
        raise SystemExit(
            f'{DATASETS_BIN!r} is not on PATH. Install with '
            "'conda install -c conda-forge ncbi-datasets-cli', "
            'or point $DATASETS at the binary.'
        ) from None
    return time.time() - t0


print('=== History-server download (small query) ===')
OUT_PATH = 'insulin_mrna.fasta'
elapsed = history_server_download(
    'nucleotide',
    # Gene-symbol field tags need the official gene symbol, not a descriptive word
    # (checked live 2026-09-19: 'insulin[Gene Name]' returns Count=0; 'INS[GENE]' does not).
    'Homo sapiens[ORGN] AND INS[GENE] AND srcdb_refseq[PROP] AND biomol_mrna[PROP]',
    OUT_PATH,
)
print(f'  elapsed: {elapsed:.1f}s')

print('\n=== Verify integrity ===')
if elapsed == 0:
    print('  No records found -- nothing to verify')
else:
    records = list(SeqIO.parse(OUT_PATH, 'fasta'))
    print(f'  {len(records)} records in file')
    for r in records[:5]:
        print(f'    {r.id}: {len(r.seq)} nt')

print('\n=== When to defect to Datasets CLI ===')
print('For genome assemblies, use:')
print("  datasets download genome accession GCF_000001405.40 --include genome,protein,gff3")
print('For all reference genomes for a taxon:')
print("  datasets download genome taxon 'Escherichia coli' --reference")
print('The Datasets CLI handles checksums, parallel download, and accession resolution.')
print('Equivalent E-utils workflow is slower and needs custom md5 handling.')

# Opt-in: a real download of a small genome, so the Datasets path is exercised
# rather than only described. Large pulls stay behind the flag -- a taxon-wide
# download is hundreds of GB and should never fire from a demo import.
if '--datasets' in sys.argv:
    print('\n  Downloading SARS-CoV-2 reference (~20 KB)...')
    elapsed = datasets_download('GCF_009858895.2', 'sars2.zip')
    print(f'  downloaded in {elapsed:.1f}s -> sars2.zip')
    print('  Inspect with: unzip -l sars2.zip')
else:
    print('\n  (pass --datasets to run a small real download)')
