'''ELink calls with explicit linkname choices; demonstrates curated vs all-hits, asymmetric round-trip, scored neighbors.'''
# Reference: biopython 1.83+, entrez direct 21.0+ | Verify API if version differs
from Bio import Entrez
import time

Entrez.email = 'your.email@example.com'
DELAY = 0.34


def linked_ids(dbfrom, db, source_id, linkname=None, cmd='neighbor'):
    kwargs = {'dbfrom': dbfrom, 'db': db, 'id': source_id, 'cmd': cmd}
    if linkname:
        kwargs['linkname'] = linkname
    h = Entrez.elink(**kwargs)
    r = Entrez.read(h); h.close()
    if not r[0]['LinkSetDb']:
        return []
    return [link['Id'] for link in r[0]['LinkSetDb'][0]['Link']]


GENE_BRCA1 = '672'

print('=== Curated vs all-protein link: ratio matters ===')
refseq = linked_ids('gene', 'protein', GENE_BRCA1, linkname='gene_protein_refseq')
time.sleep(DELAY)
all_proteins = linked_ids('gene', 'protein', GENE_BRCA1, linkname='gene_protein')
print(f'gene_protein_refseq:  {len(refseq):>5} proteins (RefSeq isoforms)')
print(f'gene_protein (all):   {len(all_proteins):>5} proteins (incl. predictions)')
print(f'Ratio: {len(all_proteins) / max(len(refseq), 1):.1f}x')
time.sleep(DELAY)

print('\n=== Asymmetric/curation-level warning: Gene -> PubMed ===')
# Note: an earlier version of this demo used PMID 35412348's pubmed_gene/pubmed_gene_rif
# links, but that PMID currently has zero links in either linkname (live NCBI, confirmed
# via cmd='acheck' -- no pubmed_gene* linkname is even offered for it). TP53's own
# gene_pubmed_rif vs gene_pubmed links reliably demonstrate the same curated-vs-all
# disparity on live data.
GENE_TP53 = '7157'
genes_via_rif = linked_ids('gene', 'pubmed', GENE_TP53, linkname='gene_pubmed_rif')
time.sleep(DELAY)
genes_via_all = linked_ids('gene', 'pubmed', GENE_TP53, linkname='gene_pubmed')
print(f'gene_pubmed_rif (curated only):      {len(genes_via_rif)} PubMed records')
print(f'gene_pubmed (text-mined + curated):  {len(genes_via_all)} PubMed records')
print(f'RIF set is a subset of the all set: {set(genes_via_rif).issubset(set(genes_via_all))}')
time.sleep(DELAY)

print('\n=== neighbor_score: relevance scores for related papers ===')
h = Entrez.elink(dbfrom='pubmed', db='pubmed', id='35412348',
                 linkname='pubmed_pubmed', cmd='neighbor_score')
r = Entrez.read(h); h.close()
if r[0]['LinkSetDb']:
    scored = [(l['Id'], int(l['Score'])) for l in r[0]['LinkSetDb'][0]['Link'][:5]]
    print(f'Top 5 by relevance score:')
    for pmid, score in scored:
        print(f'  PMID {pmid}: {score}')

print('\n=== BioProject -> SRA cross-reference ===')
h = Entrez.esearch(db='bioproject', term='PRJNA661299[BioProject]')
r = Entrez.read(h); h.close()
time.sleep(DELAY)
if r['IdList']:
    bp_uid = r['IdList'][0]
    sra_uids = linked_ids('bioproject', 'sra', bp_uid)
    print(f'PRJNA661299 -> {len(sra_uids)} SRA runs (first 5 UIDs: {sra_uids[:5]})')
