# Errors



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
- **Fix:** Use `query_raw()` (Code patterns), which raises `BioMartOutageError` for an HTML/status
  response and `BioMartResponseError` for an empty or non-TSV response. Retry only the outage error
  with backoff; inspect the query/filter values for the other error.

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
| `BioMartInputError: filter '...' is empty` | An empty list was supplied | Correct the input; no network request was made |
| `BioMartOutageError` | HTML/status page served instead of TSV | Retry with backoff; this is a service issue |
| `BioMartResponseError` | Empty, rejected, or non-TSV response | Inspect filter values and response; do not assume it is an outage |
| Timeout on big query | No filter, too many rows | Chunk by chromosome |
| Drift between re-runs | No version pinning | `useEnsembl(version=110)` |
| Row count > expected | Many-to-many cross-ref joins | Filter to canonical isoform |
| Symbol filter returns nothing | HGNC rename | Filter by Ensembl ID or HGNC ID |
| Slow on ortholog wide-table | Multi-species join expensive | Chunk by chromosome |
| `Query ERROR ... Too many attributes selected for External References` | >3 attributes from that attribute page in one query | Split into two queries, join client-side |
