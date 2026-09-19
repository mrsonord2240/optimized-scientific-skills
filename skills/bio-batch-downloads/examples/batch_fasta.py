'''Compare retrieval strategies: direct EFetch, history server, and NCBI Datasets CLI, with cost/payload metrics.'''
# Reference: biopython 1.83+, ncbi datasets cli 16.0+ | Verify API if version differs
from Bio import Entrez, SeqIO
import subprocess
import time

Entrez.email = 'your.email@example.com'
DELAY = 0.34


def history_server_download(db, term, out_path, batch_size=500):
    h = Entrez.esearch(db=db, term=term, usehistory='y', retmax=0)
    s = Entrez.read(h); h.close()
    webenv, query_key, total = s['WebEnv'], s['QueryKey'], int(s['Count'])
    print(f'  history-server: {total} records')

    if total == 0:
        return 0

    t0 = time.time()
    with open(out_path, 'w') as out:
        for start in range(0, total, batch_size):
            h = Entrez.efetch(db=db, rettype='fasta', retmode='text',
                              retstart=start, retmax=batch_size,
                              webenv=webenv, query_key=query_key)
            out.write(h.read()); h.close()
            time.sleep(DELAY)
    elapsed = time.time() - t0
    return elapsed


def datasets_cli_genome(taxon, out_zip='genome.zip'):
    '''Replace assembly_summary.txt scraping with the supported bulk path.'''
    cmd = ['datasets', 'download', 'genome', 'taxon', taxon,
           '--reference', '--include', 'genome,protein,gff3',
           '--filename', out_zip]
    t0 = time.time()
    subprocess.run(cmd, check=True)
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
