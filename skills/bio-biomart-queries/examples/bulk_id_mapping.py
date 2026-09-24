'''Bulk ID mapping via pybiomart: Ensembl Gene -> HGNC + RefSeq + UniProt.

Uses query_raw() (see SKILL.md "Querying with ID-list filters") instead of ds.query()
because pybiomart 0.2.0's Dataset.filters never recurses into id_list filter
collections, so ds.query(filters={'ensembl_gene_id': [...]}) raises BiomartException
even though ensembl_gene_id is a real, valid filter (confirmed server-side).

Limited to 3 "External References" attributes (hgnc_id, refseq_mrna, uniprotswissprot)
-- adding a 4th (e.g. entrezgene_id) makes Ensembl reject the query server-side with
"Too many attributes selected for External References". Query entrezgene_id in a
separate call and join client-side on ensembl_gene_id if you need it too.
'''
# Reference: pybiomart 0.2.0, Ensembl release 116 | checked live 2026-09-24
import sys
from pathlib import Path
from pybiomart import Server
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from biomart_query import query_raw


server = Server(host='http://www.ensembl.org')
mart = server['ENSEMBL_MART_ENSEMBL']
ds = mart['hsapiens_gene_ensembl']


ensembl_ids = [
    'ENSG00000139618',  # BRCA2
    'ENSG00000141510',  # TP53
    'ENSG00000171862',  # PTEN
    'ENSG00000146648',  # EGFR
    'ENSG00000136997',  # MYC
]

print('=== Bulk ID mapping: Ensembl Gene -> HGNC + RefSeq + UniProt ===')
df = query_raw(ds,
    attributes=[
        'ensembl_gene_id',
        'external_gene_name',
        'hgnc_id',
        'refseq_mrna',
        'uniprotswissprot',
    ],
    filters={'ensembl_gene_id': ensembl_ids},
)
print(f'  Rows: {len(df)} (note: many-to-many cross-ref joins multiply rows)')
print(df.head(10).to_string(index=False))

print('\n=== Collapse to one row per gene (most common downstream pattern) ===')
# NaN cells are floats, not None -- filter(None, x) does not drop them (NaN is
# truthy), so join on non-null values explicitly instead.
join_non_null = lambda x: ';'.join(sorted(set(str(v) for v in x if pd.notna(v))))
collapsed = (df.groupby('Gene stable ID')
               .agg({'Gene name': 'first',
                     'HGNC ID': 'first',
                     'RefSeq mRNA ID': join_non_null,
                     'UniProtKB/Swiss-Prot ID': join_non_null })
               .reset_index())
print(collapsed.to_string(index=False))


print('\n=== Discover attribute names if column labels differ ===')
print('First 10 of', len(ds.attributes), 'attributes:')
for a, info in list(ds.attributes.items())[:10]:
    print(f'  {a:<35} {info.display_name}')
