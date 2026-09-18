---
name: bio-biomart-queries
description: Bulk-query Ensembl BioMart (and other BioMart instances) for cross-database ID mapping, gene/transcript/exon coordinates, and ortholog tables. Use when batch-converting Ensembl IDs to other namespaces (HGNC, RefSeq, UniProt, Entrez), pulling gene coordinate tables for thousands of genes, building ortholog wide-tables across species, or replacing slow Ensembl REST loops with one-shot bulk export. Encodes BioMart's XML query format, R biomaRt vs Python pybiomart trade-off, mart-vs-dataset hierarchy, and the URL endpoint that's BioMart-specific (separate from rest.ensembl.org).
tool_type: mixed
primary_tool: pybiomart
license: MIT
---

## Version Compatibility

Reference examples checked live 2026-09-17 on **pybiomart 0.2.0** (the only version ever published to PyPI -- there is no 0.9 release, on any date) and R biomaRt 2.62.1 (Bioconductor 3.20); Ensembl BioMart release 116.

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show pybiomart`
- R: `packageVersion('biomaRt')`

**pybiomart 0.2.0's `Dataset.filters` only enumerates the 45 top-level filter names -- it never recurses into `id_list`-type filter collections.** `ensembl_gene_id`, `external_gene_name`, `entrezgene_id`, `hgnc_id` and similar ID-list filters are real, valid, server-side filters (confirmed by sending them directly in the martservice XML), but `ds.query(filters={'ensembl_gene_id': [...]})` raises `BiomartException: Unknown filter ensembl_gene_id` because that name never appears in `ds.filters`. This is a client-side gap in pybiomart, not an Ensembl-side removal, and discovering filters first (below) does not protect you from it -- the broken names never show up in discovery output either. Every pattern below that filters on an ID list uses the `query_raw()` workaround instead of `ds.query()` for this reason; see "Querying with ID-list filters" in Code patterns.

The BioMart XML query format is stable across Ensembl releases; the underlying mart names and attribute IDs can change between Ensembl releases. For published work, pin the Ensembl release via `useEnsembl(version=110)`.

# BioMart Queries

**"Bulk-convert IDs / pull coordinate tables / extract ortholog wide tables"** -> BioMart is the right answer for any Ensembl-rooted query producing >5,000 rows. It is a separate service from the Ensembl REST API, with separate rate behavior and a different query model (XML-based, batch-oriented). For one-off lookups (<100 records), Ensembl REST is more convenient; for bulk anything, BioMart wins.

The single most important fact: **BioMart returns a flat table from a single query**. There is no per-record loop, no rate-limit cascade, no async polling. One XML query in; one TSV out.

- Python: `pybiomart` (https://github.com/jrderuiter/pybiomart) is the lightest client
- R: `biomaRt` Bioconductor (Durinck et al. 2009 *Nat Protoc* 4:1184) is the canonical client -- more mature and Bioconductor-supported than pybiomart; prefer it for R-based pipelines
- CLI: `curl` against the XML endpoint works but is rarely used directly
- Web: `https://www.ensembl.org/biomart/martview` for interactive query design
- Non-vertebrate species: swap the host for the Ensembl Genomes BioMart, e.g. `Server(host='http://plants.ensembl.org')`

## Installation

```bash
pip install pybiomart pandas
# R:
# BiocManager::install('biomaRt')
```

## BioMart hierarchy

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
| Bulk ID mapping (>5000 IDs) | yes (1 query) | rate-limited cascade |
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

`ds.query(filters={'ensembl_gene_id': [...]})` fails as described above. The fix confirmed live
(TP53/BRCA1/PTEN/EGFR/MYC, real HGNC/RefSeq/UniProt cross-refs, checked 2026-09-17): build the same
XML `ds.query()` builds internally and send it through `ds.get()`, which skips the broken
attribute/filter-dict validation and lets Ensembl answer directly. This also guards against
Ensembl's intermittent "Service unavailable" page, which comes back as HTTP 200 and would otherwise
be silently parsed as data (see Failure modes).

```python
from io import StringIO
from xml.etree import ElementTree
import pandas as pd

def query_raw(ds, attributes, filters):
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
```

The patterns below all use `query_raw()`, defined once here, instead of `ds.query()`.

### Bulk ID mapping: Ensembl Gene -> HGNC + RefSeq + UniProt

**Goal:** Convert 5,000 Ensembl Gene IDs to HGNC symbols, RefSeq mRNA accessions, and UniProt accessions in one query.

**Approach:** `query_raw()` with three cross-ref attributes; ID list as a filter; returns one TSV.
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

ensembl_ids = ['ENSG00000139618', 'ENSG00000141510', 'ENSG00000171862']  # ...up to 5K+

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

## Failure modes

### Trying to pull >100K rows in one query
- **Trigger:** Query without any filter (e.g. all attributes for the whole human genome).
- **Mechanism:** BioMart times out or truncates on very large queries.
- **Symptom:** Empty or partial result.
- **Fix:** Chunk by chromosome; combine results client-side.

### No version pinning
- **Trigger:** `useMart('ensembl', ...)` without `version=`.
- **Mechanism:** Defaults to current release; gene model versions change quarterly.
- **Symptom:** Re-running a year later produces different rows.
- **Fix:** Pin with `useEnsembl(version=110)` or archive host URL.

### Multiple cross-refs balloon row count
- **Trigger:** Query for `ensembl_gene_id, refseq_mrna`; a gene with 10 RefSeq mRNAs produces 10 rows.
- **Mechanism:** BioMart joins on cross-refs; many-to-many produces row multiplication.
- **Symptom:** "Why do I have 50K rows for 5K input IDs?"
- **Fix:** Filter to one isoform per gene downstream; or use `ensembl_canonical` filter where available.

### Symbol-based filter misses HGNC renames
- **Trigger:** `filters={'external_gene_name': ['MARCH1']}` post-2020.
- **Mechanism:** HGNC renamed to MARCHF1; BioMart mirrors the new symbol.
- **Symptom:** Empty result for that gene.
- **Fix:** Filter by `ensembl_gene_id` or `hgnc_id`; these are stable.

### Multi-species mart query slow
- **Trigger:** Querying `mmusculus_homolog_ensembl_gene` for 30K human genes.
- **Mechanism:** Ortholog attributes are heavy; large queries take minutes.
- **Symptom:** Timeout or slow.
- **Fix:** Chunk by chromosome; or use Ensembl Compara REST for targeted lookups.

### REST loops where BioMart belongs
- **Trigger:** Loop of 5,000 Ensembl REST `/lookup/symbol` calls.
- **Mechanism:** Rate-limit cascade; 5,000 * 0.07s = 6 minutes just for the rate gate, plus HTTP overhead.
- **Symptom:** Slow; 429 errors.
- **Fix:** Switch to one BioMart query.

### Wrong mart for the question
- **Trigger:** Querying gene info from `ENSEMBL_MART_SNP`.
- **Mechanism:** SNP mart has variant attributes, not gene attributes.
- **Symptom:** Empty result or wrong fields.
- **Fix:** Discover marts with `server.marts`; pick `ENSEMBL_MART_ENSEMBL` for genes.

### Outage page silently parsed as data
- **Trigger:** Querying during an Ensembl outage window (observed live, 2026-09-17).
- **Mechanism:** Ensembl serves its `status.ensembl.org` "Service unavailable" HTML page with HTTP
  200, not an error code. `ds.query()`'s own TSV parser accepts it as valid data with no check, then
  downstream column lookups (e.g. `next(c for c in df.columns if ...)`) crash with an opaque
  `StopIteration` that gives no hint of the real cause.
- **Symptom:** `StopIteration`, or a one-row/one-column garbage DataFrame, with no BioMart error text.
- **Fix:** Use `query_raw()` (Code patterns), which checks the raw response body for an HTML/empty
  payload before parsing and raises a clear `RuntimeError` telling you to retry. Retry with backoff;
  this is a live-service availability issue, not a code error.

### ID-list filter rejected even though it's a real filter
- **Trigger:** `ds.query(filters={'ensembl_gene_id': [...]})` or `external_gene_name` / `entrezgene_id`.
- **Mechanism:** pybiomart 0.2.0's `Dataset.filters` never recurses into `id_list`-type filter
  collections (see Version Compatibility); the name is valid server-side but invisible to the client.
- **Symptom:** `BiomartException: Unknown filter ensembl_gene_id, check dataset filters for a list of
  valid filters` -- raised before any network call.
- **Fix:** Use `query_raw()` instead of `ds.query()` for any ID-list filter.

## Common errors

| Error / symptom | Cause | Solution |
|---|---|---|
| Empty result | Wrong attribute / filter name | List with `ds.attributes` and `ds.filters` |
| `BiomartException: Unknown filter ensembl_gene_id` (or `external_gene_name`, `entrezgene_id`) | ID-list filter not exposed by `Dataset.filters` in pybiomart 0.2.0 | Use `query_raw()` instead of `ds.query()` |
| `StopIteration` with no BioMart error text | Outage page served as HTTP 200, parsed as data | Use `query_raw()`; retry with backoff |
| Timeout on big query | No filter, too many rows | Chunk by chromosome |
| Drift between re-runs | No version pinning | `useEnsembl(version=110)` |
| Row count > expected | Many-to-many cross-ref joins | Filter to canonical isoform |
| Symbol filter returns nothing | HGNC rename | Filter by Ensembl ID or HGNC ID |
| Slow on ortholog wide-table | Multi-species join expensive | Chunk by chromosome |
| `Query ERROR ... Too many attributes selected for External References` | >3 attributes from that attribute page in one query | Split into two queries, join client-side |

## References

- Durinck S, Spellman PT, Birney E, Huber W. (2009) Mapping identifiers for the integration of genomic datasets with the R/Bioconductor package biomaRt. *Nat Protoc* 4:1184-1191.
- Kinsella RJ, Kahari A, Haider S, et al. (2011) Ensembl BioMarts: a hub for data retrieval across taxonomic space. *Database* 2011:bar030.
- Smedley D, Haider S, Durinck S, et al. (2015) The BioMart community portal: an innovative alternative to large, centralized data repositories. *Nucleic Acids Res* 43:W589-W598.
- pybiomart documentation: https://github.com/jrderuiter/pybiomart

## Related Skills

- ensembl-rest - Per-record Ensembl queries (BioMart's complement)
- ortholog-inference - Compara ortholog calls with confidence semantics
- uniprot-access - UniProt ID mapping (preferred for UniProt-rooted lookups and obsolete-accession resolution; BioMart is preferred for Ensembl-rooted batches >5K)
- ncbi-datasets-cli - NCBI-side bulk path for genome / gene data
- entrez-search - NCBI alternative for non-Ensembl queries
