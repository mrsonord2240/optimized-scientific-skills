'''Build a multi-species ortholog wide-table from one BioMart query; filter to 1:1 orthologs.

Uses query_raw() (see SKILL.md "Querying with ID-list filters") instead of ds.query() so
that a malformed response -- Ensembl's intermittent "Service unavailable" page, served as
HTTP 200 -- is caught with a clear RuntimeError instead of silently parsed into a garbage
DataFrame that then crashes column lookup with an opaque StopIteration.
'''
# Reference: pybiomart 0.2.0, Ensembl release 116 | checked live 2026-09-24
import sys
from pathlib import Path
from pybiomart import Server

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from biomart_query import query_raw


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
