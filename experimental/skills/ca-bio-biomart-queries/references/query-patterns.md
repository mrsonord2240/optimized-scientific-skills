# Query-Patterns



| Level | Examples |
|---|---|
| Mart | `ENSEMBL_MART_ENSEMBL` (genes), `ENSEMBL_MART_SNP` (variants), `ENSEMBL_MART_MOUSE` (mouse-specific) |
| Dataset | `hsapiens_gene_ensembl`, `mmusculus_gene_ensembl`, etc. (per species) |
| Attribute | Fields to return: `ensembl_gene_id`, `external_gene_name`, `chromosome_name`, etc. |
| Filter | Constraints on the query: `chromosome_name = 17`, `biotype = protein_coding`, etc. |

A query is: pick a mart, pick a dataset, list attributes to return, list filters to constrain. BioMart returns a single TSV.

Discovery:
```python
from pybiomart import Server
server = Server(host='http://www.ensembl.org')
print(server.marts)                                          # list marts
mart = server['ENSEMBL_MART_ENSEMBL']
print(mart.datasets)                                         # list datasets (species)
ds = mart['hsapiens_gene_ensembl']
print(ds.attributes)                                         # list attributes
print(ds.filters)                                            # list filters
```

## Decision matrix: BioMart vs Ensembl REST

| Question | BioMart | Ensembl REST |
|---|---|---|
| Bulk ID mapping (>5000 IDs) | yes (batched XML queries) | rate-limited cascade |
| Single-gene lookup | overkill | yes |
| Coordinate tables for thousands of genes | yes | rate-limited |
| Ortholog wide-table across species | yes (multi-species mart) | per-gene loop |
| VEP variant annotation | no | yes (or local VEP) |
| Sequence retrieval | partial | yes |
| Real-time | no (batch) | yes (per-record) |
| Reproducibility (version pin) | `useEnsembl(version=110)` | archive URL `e110.rest.ensembl.org` |

For >5K rows, BioMart is the right tool. For real-time per-record lookups, REST.

## Common attribute selectors

| Attribute | Returns |
|---|---|
| `ensembl_gene_id` | Stable Ensembl Gene ID |
| `ensembl_gene_id_version` | With `.N` version suffix |
| `external_gene_name` | HGNC symbol (or species-equivalent) |
| `hgnc_id`, `hgnc_symbol` | HGNC permanent ID and symbol |
| `entrezgene_id` | NCBI Gene ID |
| `refseq_mrna`, `refseq_peptide` | RefSeq accessions |
| `uniprotswissprot`, `uniprotsptrembl` | UniProt accessions |
| `chromosome_name`, `start_position`, `end_position`, `strand` | Gene coordinates |
| `transcript_count`, `exon_count` | Counts |
| `gene_biotype` | protein_coding, lncRNA, miRNA, etc. (as an *attribute*; `biotype` is the *filter* name for the same concept -- see below) |
| `description` | Free-text gene description |
| `go_id`, `name_1006`, `namespace_1003` | GO term ID, name, namespace |

## Common filter selectors

| Filter | Constraint |
|---|---|
| `ensembl_gene_id` | List of Gene IDs |
| `external_gene_name` | List of symbols |
| `entrezgene_id` | List of NCBI Gene IDs |
| `chromosome_name` | One or more chromosomes |
| `start` / `end` | Coordinate range |
| `biotype` | One or more biotypes |
| `with_<source>` | Boolean: has cross-ref to `<source>` (e.g. `with_hpa` = has Human Protein Atlas) |

## Code patterns

### Querying with ID-list filters (pybiomart 0.2.0 workaround)

`ds.query(filters={'ensembl_gene_id': [...]})` fails as described above. Use the shared
`scripts/biomart_query.py` helper: it builds the same XML `ds.query()` builds internally and sends it
through `ds.get()`, which skips the broken attribute/filter-dict validation. It also validates empty
ID lists locally, distinguishes HTML outage pages from other malformed responses, and starts one ID
list in 500-value requests before concatenating their TSV tables. If an instance responds 414
Request-URI Too Large, it halves only that failed batch automatically. Run its offline
check with `python scripts/biomart_query.py --self-test`.

```python
import sys
from pathlib import Path

sys.path.insert(0, str(Path('scripts').resolve()))
from biomart_query import query_raw
```

The patterns below all use `query_raw()` from `scripts/biomart_query.py`, instead of `ds.query()`.

### Bulk ID mapping: Ensembl Gene -> HGNC + RefSeq + UniProt

**Goal:** Convert 5,000 Ensembl Gene IDs to HGNC symbols, RefSeq mRNA accessions, and UniProt accessions in bounded BioMart batches.

**Approach:** `query_raw()` with three cross-ref attributes; ID list as a filter; returns one concatenated TSV. The helper starts at 500 IDs per request and halves a batch if that BioMart instance returns HTTP 414; it only accepts one such large list filter at a time.
Stick to 3 attributes from the "External References" attribute page (`hgnc_id`, `refseq_mrna`,
`uniprotswissprot` here) -- combining 4 or more of them (e.g. adding `entrezgene_id`) makes Ensembl
reject the query server-side with "Too many attributes selected for External References"; query
`entrezgene_id` separately and join client-side on `ensembl_gene_id` if you need it too.

**Reference (pybiomart 0.2.0, Ensembl release 116, checked 2026-09-17):**
```python
from pybiomart import Server

server = Server(host='http://www.ensembl.org')
mart = server['ENSEMBL_MART_ENSEMBL']
ds = mart['hsapiens_gene_ensembl']

ensembl_ids = ['ENSG00000139618', 'ENSG00000141510', 'ENSG00000171862']  # ...5K is split into 500-ID batches

df = query_raw(ds,
    attributes=['ensembl_gene_id', 'external_gene_name', 'hgnc_id',
                'refseq_mrna', 'uniprotswissprot'],
    filters={'ensembl_gene_id': ensembl_ids},
)
print(df.head())
# One row per (gene, cross-ref) pair; genes with multiple RefSeq mRNAs get multiple rows.
```

### Pull gene coordinate table for a chromosome

`chromosome_name` and `biotype` are top-level filters (not ID-list), so this would also work through
`ds.query()`, but `query_raw()` is used uniformly here for the same response validation. Note the
attribute is `gene_biotype` -- `biotype` is only a valid *filter* name, not an attribute name.

```python
df = query_raw(ds,
    attributes=['ensembl_gene_id', 'external_gene_name', 'chromosome_name',
                'start_position', 'end_position', 'strand', 'gene_biotype'],
    filters={'chromosome_name': '17', 'biotype': 'protein_coding'},
)
print(f'{len(df)} protein-coding genes on chr17')
```

### Bulk ortholog wide-table (human <-> mouse <-> zebrafish)

**Goal:** One TSV with human Ensembl ID, mouse ortholog Ensembl ID, zebrafish ortholog Ensembl ID per row.

**Approach:** Ortholog attributes from the human mart query both species' orthologs.

```python
df = query_raw(ds,
    attributes=['ensembl_gene_id', 'external_gene_name',
                'mmusculus_homolog_ensembl_gene', 'mmusculus_homolog_orthology_type',
                'drerio_homolog_ensembl_gene', 'drerio_homolog_orthology_type'],
    filters={'chromosome_name': '17'},
)
# pybiomart columns use the mart display names, which can vary across releases.
# Resolve column names defensively rather than hardcoding strings:
mouse_type_col = next(c for c in df.columns if 'Mouse' in c and 'type' in c)
zebra_type_col = next(c for c in df.columns if 'Zebrafish' in c and 'type' in c)
df_one2one = df[(df[mouse_type_col] == 'ortholog_one2one') &
                (df[zebra_type_col] == 'ortholog_one2one')]
print(f'{len(df_one2one)} 1:1 orthologs across all three species on chr17')
```

### GO term annotation for a gene set

```python
df = query_raw(ds,
    attributes=['ensembl_gene_id', 'external_gene_name',
                'go_id', 'name_1006', 'namespace_1003'],
    filters={'external_gene_name': ['TP53', 'BRCA1', 'MYC', 'EGFR']},
)
# Long format: one row per (gene, GO term) pair
```

### Version-pinned query (R biomaRt)

```r
# Reference: Bioconductor biomaRt 2.58+ | Verify API if version differs
library(biomaRt)

# Pin to release 110 for reproducibility
ensembl <- useEnsembl(biomart='genes', dataset='hsapiens_gene_ensembl', version=110)

# Or via host URL (for older or specific assemblies)
# ensembl <- useMart('ENSEMBL_MART_ENSEMBL',
#                     dataset='hsapiens_gene_ensembl',
#                     host='https://nov2020.archive.ensembl.org')

df <- getBM(
    attributes = c('ensembl_gene_id', 'external_gene_name', 'entrezgene_id',
                   'uniprotswissprot', 'refseq_mrna'),
    filters = 'ensembl_gene_id',
    values = c('ENSG00000139618', 'ENSG00000141510'),
    mart = ensembl
)
head(df)
```

### Discover attributes / filters programmatically

```python
# What attributes are available?
attrs = ds.attributes
ortho_attrs = [a for a in attrs if 'homolog' in a]
print(f'{len(ortho_attrs)} ortholog attributes; first 5: {ortho_attrs[:5]}')

# What filters?
filts = ds.filters
chrom_filts = [f for f in filts if 'chrom' in f]
```
