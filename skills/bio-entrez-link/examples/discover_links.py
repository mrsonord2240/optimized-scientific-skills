'''Enumerate every linkname available for a (dbfrom, source-id) pair via cmd=acheck.'''
# Reference: biopython 1.83+, entrez direct 21.0+ | Verify API if version differs
from Bio import Entrez
import time

Entrez.email = 'your.email@example.com'
DELAY = 0.34


def discover_links(dbfrom, source_id):
    '''Return (LinkName, DbTo, MenuTag) tuples for the (dbfrom, source-id) pair.

    LinkInfo entries are keyed 'LinkName' (not 'Name') on current Biopython/NCBI,
    and some entries (e.g. nuccore_nuccore_mrnaonly, pubmed_pmc_local) omit
    MenuTag/HtmlTag entirely -- both must be guarded, not indexed directly.
    '''
    h = Entrez.elink(dbfrom=dbfrom, id=source_id, cmd='acheck')
    r = Entrez.read(h); h.close()
    info = r[0]['IdCheckList']['IdLinkSet'][0]['LinkInfo']
    # Schema guard: fail loudly if NCBI/Biopython renames or drops a key (the 'Name' ->
    # 'LinkName' drift was silent until it raised KeyError deep in a pipeline).
    for i in info:
        assert 'LinkName' in i and 'DbTo' in i, f'LinkInfo schema changed, keys={sorted(i)}'
    return [(i['LinkName'], i['DbTo'], i.get('MenuTag', '<none>')) for i in info]


def main():
    print('=== Gene 672 (BRCA1): all available link tables ===')
    for name, target, label in discover_links('gene', '672'):
        print(f'  {name:<42} -> {target:<14} ({label})')
    time.sleep(DELAY)

    print('\n=== Nucleotide UID 1732746264 (NM_007294.4): link tables ===')
    for name, target, label in discover_links('nucleotide', '1732746264'):
        print(f'  {name:<42} -> {target:<14} ({label})')
    time.sleep(DELAY)

    print('\n=== PubMed 35412348: link tables ===')
    for name, target, label in discover_links('pubmed', '35412348'):
        print(f'  {name:<42} -> {target:<14} ({label})')


if __name__ == '__main__':
    main()
