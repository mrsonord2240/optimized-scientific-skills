---
name: bio-uniprot-access
category: Data Analysis
description: Query UniProt's REST API (post-2022 endpoint at rest.uniprot.org) for protein sequences, annotations, GO terms, cross-references, ID mappings, and proteomes. Use when fetching UniProtKB entries, navigating the JSON schema, choosing between UniProtKB/UniRef/UniParc/Proteomes resources, deciding stream vs search endpoint for batch retrieval, running ID-mapping jobs with the async pattern, handling isoform suffixes, or filtering reviewed Swiss-Prot vs auto-annotated TrEMBL. Encodes the legacy URL migration (2022), the new JSON schema layout, and bulk-pull patterns.
tool_type: python
primary_tool: requests
license: MIT
author: GPTomics
---

## Version Compatibility

Reference code needs requests 2.31+ and pandas 2.2+ (checked with 2.34.2 and 2.3.3). The live service reported UniProt release 2026_03 (`X-UniProt-Release` header, 2026-10-02).

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show requests pandas`
- API surface: confirm endpoint URLs match https://www.uniprot.org/help/api

The REST API JSON schema is stable within a release; major schema changes are documented at https://www.uniprot.org/release-notes. The 2022 migration broke the legacy `https://www.uniprot.org/uniprot/...` endpoints. Pin the UniProt release (from the `X-UniProt-Release` response header) in any citation.

# UniProt Access

**"Get protein information from UniProt"** -> Two facts dominate every UniProt workflow: (1) **the API endpoint migrated in 2022** from `https://www.uniprot.org/uniprot/...` to `https://rest.uniprot.org/uniprotkb/...` with a substantially different JSON schema; pre-2022 code does not work as-is. (2) **`?fields=`** is essential -- default JSON returns the full entry (~20-30 KB each); for bulk pulls, request only the fields actually needed.

The databases under the UniProt umbrella have different scopes:

- **UniProtKB**: the curated knowledgebase -- Swiss-Prot (manually reviewed, ~576K entries in release 2026_03) + TrEMBL (auto-annotated, ~149M in 2026_03). Always specify `reviewed:true` for high-quality reference work.
- **UniRef**: clustered sequences at 100%, 90%, 50% identity. UniRef50 is the standard for redundancy reduction.
- **UniParc**: archival "every unique sequence ever seen" -- for provenance and historical lookup.
- **Proteomes**: organism-level groupings; reference proteomes (one per species) are the canonical subset.

- Python: `scripts/uniprot_client.py` (wraps `requests`); `Bio.ExPASy.get_sprot_raw()` for legacy SwissProt flat files
- CLI: `curl https://rest.uniprot.org/uniprotkb/P04637.json`

No API key required. Rate limit is generous (~200 req/sec tolerated empirically); ID mapping has its own job queue. The client retries 429/5xx with `Retry-After` and uses a 60 s per-request timeout.

## Client and examples

`scripts/uniprot_client.py` is the single authoritative client (`sys.path.insert(0, 'scripts')`): `fetch_entry_json`, `search_tsv`, `stream_tsv`, `map_ids` (async, hard timeout, reads the full `/idmapping/stream` result and lists unmapped inputs in `failedIds`), `resolve_obsolete`, `list_isoforms`, `fetch_isoform_fasta`, `xref_summary`, `download_proteome`, `uniref_cluster`.

| Task | Example |
|---|---|
| Entry parse, bulk TSV search, stream, ID mapping | `examples/uniprot_query.py` |
| Isoforms, cross-reference summary, proteome download | `examples/isoforms_and_xrefs.py` |

```python
from uniprot_client import search_tsv, map_ids

df, total = search_tsv('organism_id:9606 AND reviewed:true AND keyword:KW-0418',
                       fields=['accession', 'gene_primary', 'protein_name', 'length', 'xref_pdb'],
                       with_total=True)   # total = server match count; warns if df is a truncated page
mapping = map_ids(['ENSG00000141510'], from_db='Ensembl', to_db='UniProtKB-Swiss-Prot')
for r in mapping['results']: print(r['from'], '->', r['to'])
for failed in mapping['failedIds']: print(failed, 'NOT MAPPED')
```

## Endpoint reference

Base: `https://rest.uniprot.org/`

| Resource | Endpoint | Use |
|---|---|---|
| Single entry | `/uniprotkb/{accession}` | One protein record |
| Search | `/uniprotkb/search` | Query with up to 500 results per page |
| Stream | `/uniprotkb/stream` | No 500-result limit; for bulk |
| Batch by accession | `/uniprotkb/accessions` | Multiple specific accessions |
| ID Mapping (run) | `/idmapping/run` | Submit conversion job |
| ID Mapping (status) | `/idmapping/status/{jobId}` | Poll |
| ID Mapping (results) | `/idmapping/results/{jobId}` | One page only (default 25 rows) |
| ID Mapping (stream) | `/idmapping/stream/{jobId}` | All rows; use this |
| UniRef entry | `/uniref/{cluster_id}` | One cluster |
| UniRef search | `/uniref/search` | UniRef cluster queries |
| Proteome | `/proteomes/{upid}` | Organism proteome |
| Proteome FASTA | `/uniprotkb/stream?query=proteome:{upid}&format=fasta&compressed=true` | Download whole proteome, reviewed and unreviewed TrEMBL entries alike (human UP000005640: 147,520 entries, 37.8 MB gzip); add `AND reviewed:true` to the query for Swiss-Prot only (`/proteomes/{upid}.fasta.gz` returns 400) |
| Taxonomy | `/taxonomy/{taxid}` | Taxonomy info |

Append `.json`, `.fasta`, `.tsv`, `.xml`, `.txt`, or `.gff` to single-entry URLs to control format.

## Search query syntax

UniProt search queries use a Lucene-like syntax distinct from Entrez:

| Query | Means |
|---|---|
| `gene:TP53` | Gene name TP53 |
| `gene_exact:TP53` | Exact gene name (no wildcard match) |
| `organism_id:9606` | Human (NCBI taxonomy ID) |
| `organism_name:"Homo sapiens"` | By name (slower than taxid) |
| `reviewed:true` | Swiss-Prot only |
| `reviewed:false` | TrEMBL only |
| `length:[100 TO 500]` | Sequence length range |
| `go:0006915` | GO term (apoptosis) |
| `keyword:KW-0067` | UniProt keyword; use the KW id (`keyword:"Kinase"` matches more than `keyword:KW-0418`) |
| `ec:2.7.1.1` | Enzyme classification |
| `database:pdb` | Has PDB cross-ref (`xref:pdb` silently matches 0) |
| `existence:1` | Evidence at protein level (1 = strongest) |

Combine: `organism_id:9606 AND reviewed:true AND keyword:KW-0067 AND database:pdb`.

## `?fields=` for bulk pulls

Restrict fields for batch work; the TSV form is the compact one (`format=tsv&fields=...&size=500`).

| Field | Returns |
|---|---|
| `accession`, `id` | Primary accession (P04637), entry name (P53_HUMAN) |
| `gene_names`, `gene_primary` | All gene names / primary only |
| `protein_name` | Recommended name |
| `organism_name`, `organism_id` | Species |
| `length`, `mass`, `sequence` | Sequence stats and the sequence |
| `cc_function`, `cc_subcellular_location` | Function and localization comments |
| `ft_domain`, `ft_binding`, `ft_act_site` | Domain/site features |
| `go_p`, `go_c`, `go_f` | GO biological process / cellular component / molecular function |
| `xref_pdb`, `xref_alphafolddb`, `xref_ensembl`, `xref_refseq` | Cross-references |
| `keyword`, `ec`, `reviewed` | Keywords, enzyme class, Swiss-Prot flag |
| `cc_alternative_products` | Isoforms |

Full list: https://rest.uniprot.org/configure/uniprotkb/result-fields.

## Stream vs search vs accessions

| Endpoint | When | Limit |
|---|---|---|
| `/uniprotkb/{acc}` | One accession | 1 entry |
| `/uniprotkb/accessions?accessions=...` | Several known accessions | Up to ~100 per call |
| `/uniprotkb/search?query=...` | Query-driven; need pagination | 500 per page; `cursor=` for paging |
| `/uniprotkb/stream?query=...` | Bulk query (>500) | No hard limit; one HTTP stream |

For 1000+ results use `/stream`. A search that matches 800 records and reads only the first page silently drops the tail; `X-Total-Results` carries the true count (`search_tsv(..., with_total=True)`).

## JSON schema navigation (post-2022 layout)

The schema is deeply nested; use `.get()` chains because many fields are optional.

| Want | Path |
|---|---|
| Accession / entry name | `primaryAccession` / `uniProtkbId` |
| Sequence, length | `sequence.value`, `sequence.length` |
| Recommended name | `proteinDescription.recommendedName.fullName.value` |
| Primary gene | `genes[0].geneName.value` |
| Cross-references | `uniProtKBCrossReferences[]` with `database` and `id` |
| Domains / binding sites | `features[]` filtered on `type` (`Domain`, `Binding site`) |
| Isoforms | `comments[]` where `commentType == 'ALTERNATIVE PRODUCTS'`, then `isoforms[].name.value` |
| Reviewed flag | `entryType` starts with `UniProtKB reviewed` |

## Isoform handling

The bare accession (`P04637`) returns the canonical sequence. Isoforms use `-2`, `-3` suffixes (`P04637-2`); fetch with `/uniprotkb/P04637-2.fasta`. `list_isoforms()` reads the `ALTERNATIVE PRODUCTS` comment; iterate and fetch each isoform separately.

## ID Mapping API (async)

1. **Submit**: `POST /idmapping/run` with `ids`, `from`, `to`.
2. **Poll**: `GET /idmapping/status/{jobId}` -- returns `{'jobStatus': 'RUNNING'}` until done.
3. **Fetch**: `GET /idmapping/stream/{jobId}`. `/idmapping/results/{jobId}` is paginated (25 rows by default) and silently drops the rest if read once.

Jobs usually finish in ~30 s; larger batches take 5-10 min. **Always set a poll timeout** -- the API does not fail stuck jobs, and a network glitch can leave status at RUNNING forever (`map_ids` raises `TimeoutError`, and `RuntimeError` if the job reports `FAILED`).

`to=UniProtKB` returns TrEMBL rows too (e.g. 18 rows for TP53 Ensembl); use `to=UniProtKB-Swiss-Prot` for reviewed entries only. Inputs with no mapping are returned in `failedIds` by `map_ids` -- check it.

| From | To | Notes |
|---|---|---|
| `UniProtKB_AC-ID` | `UniProtKB` | Resolve obsolete to current accessions |
| `Gene_Name` | `UniProtKB` | Symbol -> accession (lossy; check matches) |
| `Ensembl` | `UniProtKB` | Ensembl Gene/Transcript/Protein |
| `EMBL-GenBank-DDBJ` | `UniProtKB` | INSDC nucleotide accessions |
| `RefSeq_Protein` | `UniProtKB` | NP_/XP_ accessions |
| `PDB` | `UniProtKB` | PDB chain to protein |
| `UniProtKB` | `EMBL-GenBank-DDBJ` | Reverse direction |
| `Ensembl` | `UniProtKB-Swiss-Prot` | Reviewed entries only |

Full from/to list at https://rest.uniprot.org/configure/idmapping/fields.

## Failure modes

- **Legacy URL still in code** -- `https://www.uniprot.org/uniprot/{acc}.json` gives 404 or `KeyError` from old field paths. Use `rest.uniprot.org/uniprotkb/{acc}.json` and the nested layout.
- **`fields=` omitted on bulk** -- 1000 full entries are 20-30 MB; slow, memory-heavy, rate-limit prone.
- **500-record cap** -- `/search` returns one page; use `/stream` or paginate with `cursor`. `search_tsv` warns when truncated.
- **ID-mapping results truncated** -- reading `/idmapping/results` once returns 25 rows; use `/idmapping/stream` (as `map_ids` does).
- **ID-mapping poll hangs** -- always pass `timeout=`.
- **Isoforms missed** -- the default fetch returns canonical only; fetch `-N` accessions.
- **TrEMBL flood** -- searches without `reviewed:true` return millions of auto-annotated hits.
- **Obsolete accessions** -- secondary (merged) accessions redirect to the primary entry; deleted ones return HTTP 200 with `entryType: Inactive` and no sequence (`fetch_entry_json` raises `ValueError` with the reason). Resolve with `UniProtKB_AC-ID -> UniProtKB`.
- **Gene-symbol ambiguity** -- `gene:TP53` matches all species; add `organism_id:9606`, or use `gene_exact:`.

## Common errors

| Error / symptom | Cause | Solution |
|---|---|---|
| 404 on legacy URL | Pre-2022 endpoint | Use rest.uniprot.org/uniprotkb/ |
| `KeyError` on old field path | Schema migration 2022 | Update to nested layout; use `.get()` |
| Bulk fetch very slow | Default JSON entry size | Specify `fields=` for TSV bulk |
| Mid-pagination data missing | 500-record cap | Use /stream or paginate with cursor |
| ID mapping job hangs | API doesn't fail stuck jobs | Set `timeout=` on poll loop |
| Mapping returns 25 rows / an input missing | Read first results page only | Use `/idmapping/stream`; check `failedIds` |
| HTTP 400 on `fields=` | Invalid field name (e.g. `ft_active_site`) | Use names from `/configure/uniprotkb/result-fields` (`ft_act_site`) |
| Mixed-species search results | Symbol shared across species | Add `organism_id:` filter |
| Missing isoform | Default returns canonical only | Fetch with `-N` suffix per isoform |

## References

- The UniProt Consortium. (2024) UniProt: the Universal Protein Knowledgebase in 2025. *Nucleic Acids Res* 53:D609-D617.
- Bursteinas B, Britto R, Bely B, et al. (2016) Minimizing proteome redundancy in the UniProt Knowledgebase. *Database* 2016:baw139.
- UniProt help: https://www.uniprot.org/help/api
- UniProt REST: https://rest.uniprot.org

## Related Skills

- entrez-fetch - NCBI protein records (RefSeq, GenPept) alternative
- biomart-queries - Alternative ID-mapping path via BioMart (preferred for Ensembl-rooted batches >5K; UniProt /idmapping/run is preferred for obsolete-accession resolution and any UniProt-rooted mapping)
- ortholog-inference - Resolve UniProt accessions used by OMA orthology queries
- structural-biology/structure-io - Download PDB structures referenced from UniProt
- structural-biology/alphafold-predictions - AlphaFoldDB entries cross-referenced in UniProt
- pathway-analysis/go-enrichment - Use GO annotations pulled from UniProt
