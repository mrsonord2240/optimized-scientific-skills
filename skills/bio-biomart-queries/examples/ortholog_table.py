'''Build a multi-species ortholog wide-table from one BioMart query; filter to 1:1 orthologs.

Uses query_raw() (see SKILL.md "Querying with ID-list filters") instead of ds.query() so
that a malformed response -- Ensembl's intermittent "Service unavailable" page, served as
HTTP 200 -- is caught with a clear RuntimeError instead of silently parsed into a garbage
DataFrame that then crashes column lookup with an opaque StopIteration.
'''
# Reference: pybiomart 0.2.0, Ensembl release 116 | checked live 2026-09-17
from io import StringIO
from xml.etree import ElementTree
from pybiomart import Server
import pandas as pd


def query_raw(ds, attributes, filters):
    '''Bypass Dataset.query()'s filter-name validation via ds.get() and validate
    the raw response before parsing.'''
    root = ElementTree.Element('Query')
    root.set('virtualSchemaName', 'default')
    root.set('formatter', 'TSV')
    root.set('header', '1')
    root.set('uniqueRows', '1')
    root.set('datasetConfigVersion', '0.6')
    dataset_el = ElementTree.SubElement(root, 'Dataset')
    dataset_el.set('name', ds.name)
    dataset_el.set('interface', 'default')
    for name, value in filters.items():
        f = ElementTree.SubElement(dataset_el, 'Filter')
        f.set('name', name)
        f.set('value', ','.join(value) if isinstance(value, (list, tuple)) else str(value))
    for name in attributes:
        a = ElementTree.SubElement(dataset_el, 'Attribute')
        a.set('name', name)

    response = ds.get(query=ElementTree.tostring(root))
    body = response.text.strip()
    if 'Query ERROR' in body:
        raise RuntimeError(f'BioMart rejected the query: {body}')
    if body.lower().startswith('<html') or not body:
        raise RuntimeError(
            'BioMart returned a non-TSV response (an outage page served with HTTP '
            '200, or an empty body) -- retry with backoff, this is not a code error.'
        )
    return pd.read_csv(StringIO(body), sep='\t')


server = Server(host='http://www.ensembl.org')
mart = server['ENSEMBL_MART_ENSEMBL']
ds = mart['hsapiens_gene_ensembl']


print('=== Ortholog wide-table: human + mouse + zebrafish (chr17 subset) ===')
df = query_raw(ds,
    attributes=[
        'ensembl_gene_id',
        'external_gene_name',
        'mmusculus_homolog_ensembl_gene',
        'mmusculus_homolog_orthology_type',
        'drerio_homolog_ensembl_gene',
        'drerio_homolog_orthology_type',
    ],
    filters={'chromosome_name': '17'},
)
print(f'  All chr17 genes with any ortholog row: {len(df)}')
print(df.head(8).to_string(index=False))


print('\n=== 1:1 across all three species ===')
# Column labels are the BioMart display names; check df.columns to confirm.
mouse_col = next(c for c in df.columns if 'Mouse' in c and 'type' in c)
zebra_col = next(c for c in df.columns if 'Zebrafish' in c and 'type' in c)
one2one = df[(df[mouse_col] == 'ortholog_one2one') &
             (df[zebra_col] == 'ortholog_one2one')]
print(f'  1:1 in mouse AND zebrafish: {len(one2one)}')
print(one2one.head(10).to_string(index=False))


print('\n=== Compare to per-gene Ensembl REST cost ===')
print(f'  {len(df)} chr17 genes via BioMart: one query, no rate limit')
print(f'  Same via Ensembl REST /homology/symbol: {len(df)} calls * 0.07s sleep = '
      f'{len(df) * 0.07 / 60:.1f} min minimum + HTTP overhead')
print(f'  Use REST only for <100 genes or real-time queries.')
