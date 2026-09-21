---
name: bio-entrez-search
description: Search NCBI databases using Biopython Bio.Entrez (ESearch, EInfo, ESpell), including cross-database counts via an ESearch loop (EGQuery is broken on current Biopython/NCBI — see below). Use when finding records by keyword, building reproducible field-qualified queries, navigating the Entrez Query Translator, exploiting the history server for large result sets, handling retmax caps, or interpreting weekly index lag. Covers PubMed, Nucleotide, Protein, Gene, SRA, GEO, Assembly, Taxonomy, ClinVar, dbSNP.
tool_type: python
primary_tool: Bio.Entrez
license: MIT
author: GPTomics
---

## Version Compatibility

Reference examples tested with: BioPython 1.83+, Entrez Direct 21.0+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show biopython` then `help(Bio.Entrez.esearch)` to check signatures
- CLI: `esearch -version` then `esearch -help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Entrez Search

**"Find NCBI records matching a query"** -> ESearch returns matching record UIDs (not full records) from one NCBI database; looping ESearch over a curated database list returns counts across all databases (`Bio.Entrez.egquery()` is broken — see "Cross-database counts" below); EInfo describes a database's searchable fields and update timestamp.

The single most important fact: ESearch returns *UIDs* (PMIDs, GI numbers, gene IDs, etc.), not records. To get content the agent must call EFetch or ESummary. Forgetting this is the most common Entrez mistake.

- Python: `Entrez.esearch(db=..., term=...)` (BioPython)
- CLI: `esearch -db pubmed -query 'CRISPR[Title]'` (Entrez Direct, NBK179288)
- R: `entrez_search(db=..., term=...)` (rentrez)

## Required Setup

```python
from Bio import Entrez
import os
import time

Entrez.email = 'researcher@institution.edu'       # NCBI requires; sets User-Agent
Entrez.api_key = os.environ.get('NCBI_API_KEY')   # 3 -> 10 req/sec; get at ncbi.nlm.nih.gov/account/settings/
Entrez.tool = 'project-name'                      # appears in NCBI usage logs; helps if rate-throttled
```

Never hardcode a real key in source. `Entrez.api_key` is `None` if the env var is unset, which Biopython
treats the same as not passing a key (falls back to the 3 req/sec ceiling).

**`term` strings are URL query parameters, not code** — Biopython URL-encodes them automatically (do
not pre-encode); there is no shell/eval injection risk. Do sanity-check length and count before a
large batch loop: an extremely long term (hundreds of IDs `OR`-joined into one string) is a common,
silent way to trip `HTTPError 400` — see "Common errors" below. Chunk large ID sets through EPost
(200 IDs/call, see "retmax silent caps") instead of building one giant `term`.

## What ESearch actually does

ESearch sends the query string through the **Entrez Query Translator (EQT)**, which rewrites unqualified terms into the canonical `term[field]` form, then runs the rewritten query against the per-database index. The result is a list of UIDs plus a `QueryTranslation` string showing exactly what was searched. Reproducible work always inspects `QueryTranslation` and builds queries that are translation-stable from the start.

```python
handle = Entrez.esearch(db='nucleotide', term='human BRCA1')
record = Entrez.read(handle)
handle.close()
print(record['QueryTranslation'])
# '("homo sapiens"[Organism] OR human[All Fields]) AND (BRCA1[Gene Name] OR BRCA1[All Fields])'
```

The translator may expand `human` to the full taxonomy subtree, or coerce a gene symbol to `[All Fields]` if the symbol isn't unambiguous. Use field-qualified terms (`Homo sapiens[ORGN] AND BRCA1[Gene Name]`) for any query that will be re-run later.

## Decision table: which utility for which question

| Question | Utility | Returns | Cost |
|---|---|---|---|
| "How many records match X in PubMed?" | ESearch with `retmax=0` | Count + WebEnv | 1 call |
| "Give me 20 matching UIDs" | ESearch | UIDs | 1 call |
| "Give me ALL matching UIDs (>10K)" | ESearch + `usehistory='y'` | WebEnv/QueryKey | 1 call (then EFetch chunks server-side) |
| "Does record X exist in db Y?" | ESearch with `term='X[Accn]'` | UIDs | 1 call |
| "Which of the 10 `CURATED_DBS` databases mention X?" | ESearch loop over `CURATED_DBS` (`retmax=0`) — see `examples/global_query.py` | Counts per db | N calls, ~N x 0.34s |
| "What searchable fields does db Y have?" | EInfo with `db=Y` | FieldList | 1 call |
| "Last update timestamp for db Y?" | EInfo with `db=Y` | `LastUpdate` | 1 call |
| "Did the user misspell X?" | ESpell | Spelling suggestion | 1 call |

### Cross-database counts (EGQuery is broken — use the ESearch loop)

`Bio.Entrez.egquery()` does not exist on Biopython >=1.85: `AttributeError: module 'Bio.Entrez' has
no attribute 'egquery'` (confirmed on the installed 1.88). `Bio/Entrez/__init__.py`'s own module
docstring still lists `egquery` as a provided utility, but the function was dropped from the public
API without a docstring update. Calling the underlying NCBI endpoint directly is not a working
workaround either: as of 2026-09-17 `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/egquery.fcgi`
301-redirects to `ext-http-eutils.linkerd.ncbi.nlm.nih.gov`, an internal-only NCBI hostname that does
not resolve outside their network.

Use the documented fallback instead: loop ESearch with `retmax=0` over a curated database list. It's
slower (N calls instead of 1) but never lags the per-database index the way EGQuery's counts could.

**`CURATED_DBS` covers 10 of the 38 databases EInfo currently lists (confirmed live 2026-09-17) — it
is not exhaustive.** EGQuery, which it replaces, really did query every database in one call; this
loop only checks the 10 named below. To widen it, pass any EInfo-listed db name(s) in — e.g.
`cross_db_counts(term, dbs=Entrez.read(Entrez.einfo())['DbList'])` — at the cost of one ESearch call
per database instead of 10.

See `examples/global_query.py` for the runnable pattern, reproduced here:

```python
CURATED_DBS = ['pubmed', 'pmc', 'nucleotide', 'protein', 'gene', 'sra', 'gds', 'bioproject', 'biosample', 'clinvar']

def cross_db_counts(term, dbs=CURATED_DBS):
    counts = {}
    for db in dbs:
        h = Entrez.esearch(db=db, term=term, retmax=0)
        r = Entrez.read(h); h.close()
        counts[db] = int(r['Count'])
        time.sleep(0.34)
    return counts
```

## retmax silent caps

| Endpoint behavior | Cap | Workaround |
|---|---|---|
| Default `retmax` | 20 | Set explicitly |
| Legacy esearch.fcgi (no `usehistory`) | **9,999** silent cap | Use history server |
| `usehistory='y'` + ESearch | 100,000 per page | Page with `retstart` against the WebEnv |
| EPost (to push IDs server-side) | 200 IDs per call | Chunk to multiple EPost calls; union with QueryKey |

The 9,999 cap is the bug that has shipped in countless lab pipelines: query returns "Count: 78,432" but `IdList` has 9,999 entries and there is no error. Always set `retmax` explicitly and either page or move to `usehistory='y'` whenever `Count > retmax`.

## History server (WebEnv/QueryKey) semantics

| Property | Value |
|---|---|
| TTL | 8 hours absolute (per NCBI E-utils help, 2024) |
| Idle eviction | Empirically ~15 min under load; can be shorter |
| Chaining | Run another ESearch against `WebEnv` with `term='#1 AND #2'` to intersect prior QueryKeys |
| Persistence | Session is per WebEnv string; do NOT share across processes when isolation matters |
| Failure mode | Expired session returns HTTP 200 with `<ERROR>WebEnv not found</ERROR>` — must parse body, not status |

Chaining example:
```python
h1 = Entrez.esearch(db='pubmed', term='CRISPR[Title]', usehistory='y')
r1 = Entrez.read(h1); h1.close()
webenv = r1['WebEnv']

h2 = Entrez.esearch(db='pubmed', term='2024[PDAT]', usehistory='y', WebEnv=webenv)
r2 = Entrez.read(h2); h2.close()

# Intersect QueryKey #1 (CRISPR) AND #2 (2024) into a new key
h3 = Entrez.esearch(db='pubmed', term=f'#{r1["QueryKey"]} AND #{r2["QueryKey"]}',
                    usehistory='y', WebEnv=webenv)
r3 = Entrez.read(h3); h3.close()
print(f'CRISPR & 2024: {r3["Count"]}')
```

## Index lag

NCBI's Entrez indexer runs nightly (US Eastern). Records submitted Monday morning typically appear in ESearch results Wednesday at earliest. PubMed has additional MEDLINE indexing lag (1-3 weeks for full MeSH terms). For freshly-deposited data the more reliable check is EFetch on the known accession or NCBI Datasets API for genomes.

## Field-qualified query patterns (per database)

| Database | Common fields | Notes |
|---|---|---|
| pubmed | `[Title]`, `[TIAB]` (title+abstract), `[MeSH]`, `[Author]`, `[Journal]`, `[PDAT]`, `[DCOM]`, `[PMC]` | `[TIAB]` is more permissive than `[Title]`; `[MeSH]` requires the term to be indexed (lags); PMC full-text subset is `pubmed pmc[sb]`, not a separate db — the underlying UIDs are still PubMed's |
| nucleotide | `[Organism]`, `[Gene Name]`, `[Accn]`, `[SLEN]`, `[Filter]`, `[PROP]` | `srcdb_refseq[PROP]` restricts to RefSeq; `biomol_genomic[PROP]` filters molecule type |
| protein | `[Organism]`, `[Gene Name]`, `[Accn]`, `[MOLWT]`, `[PROP]` | `swissprot[Filter]` restricts to reviewed |
| gene | `[Gene/Locus]`, `[Organism]`, `[Chromosome]`, `[Gene Type]` | `[Gene Type]` includes `protein-coding`, `pseudo`, `ncRNA` |
| sra | `[Organism]`, `[Platform]`, `[Strategy]`, `[Library Source]`, `[BioProject]` | `[Strategy]` accepts `RNA-Seq`, `WGS`, `ChIP-Seq`, etc. |
| gds (GEO) | `[Organism]`, `[Entry Type]`, `[GDS Type]`, `[Platform]` | `gse[Entry Type]` for Series, `gds[Entry Type]` for curated DataSets |
| taxonomy | `[Scientific Name]`, `[Common Name]`, `[Rank]`, `[TXID]` | TXID is the numeric taxonomy ID |
| clinvar | `[Gene Name]`, `[Clinical Significance]`, `[Variation Type]` | `pathogenic[CLIN]` for pathogenic only |

### Filter properties that newcomers miss

```python
# Curated RefSeq mRNA only, human, between 500 and 5000 nt
term = 'Homo sapiens[ORGN] AND srcdb_refseq[PROP] AND biomol_mrna[PROP] AND 500:5000[SLEN]'

# Reviewed SwissProt human kinases
term = 'Homo sapiens[ORGN] AND swissprot[Filter] AND kinase[Protein Name]'

# PubMed: human studies in last 30 days, full-text in PMC
term = 'CRISPR[Title] AND humans[MeSH Terms] AND last 30 days[EDAT] AND pubmed pmc[sb]'
```

### Organism field gotcha

`[Organism]` (and the alias `[ORGN]`) is **taxonomy-walked**: searching `mammalia[ORGN]` returns records from every species in Mammalia. To get records tagged at exactly that node use `[Organism:exp]` (no taxonomic expansion). Most workflows want the default walk, but multi-species queries that "blow up" by 100x are almost always a missing `:exp`.

## Code patterns

Runnable versions of the ESearch/EInfo/cross-database patterns below live in `examples/basic_search.py`,
`examples/database_info.py`, and `examples/global_query.py` respectively.

### Single search with explicit retmax

**Goal:** Get matching UIDs for a focused query without hitting silent caps.

**Approach:** Set `retmax` explicitly to the maximum the caller wants; if `Count > retmax` either page or switch to history server.

**Reference (BioPython 1.83+):**
```python
def search_ncbi(db, term, max_results=100):
    handle = Entrez.esearch(db=db, term=term, retmax=max_results)
    record = Entrez.read(handle); handle.close()
    count = int(record['Count'])
    if count > max_results:
        print(f'WARNING: {count} matched, returning first {max_results}; use history server for full set')
    return record['IdList'], count, record['QueryTranslation']
```

### Paged retrieval (only when enumeration without fetching is required)

**Goal:** Stream all matching UIDs to a file when downstream work can't use the history server.

**Approach:** Page through `retstart` increments; respect rate limit; stop at total.

```python
def stream_all_ids(db, term, batch_size=10000):
    h = Entrez.esearch(db=db, term=term, retmax=0)
    total = int(Entrez.read(h)['Count']); h.close()
    delay = 0.1 if Entrez.api_key else 0.34
    for start in range(0, total, batch_size):
        h = Entrez.esearch(db=db, term=term, retstart=start, retmax=batch_size)
        r = Entrez.read(h); h.close()
        for uid in r['IdList']:
            yield uid
        time.sleep(delay)
```

For any download workflow, history-server retrieval is strictly better — see `batch-downloads` skill.

### History server for downstream EFetch

**Goal:** Push a large result set to NCBI servers so EFetch can pull it in batches without re-sending IDs.

**Approach:** ESearch with `usehistory='y'`; capture WebEnv and QueryKey; pass to EFetch.

**Reference (BioPython 1.83+):**
```python
h = Entrez.esearch(db='nucleotide',
                   term='Homo sapiens[ORGN] AND srcdb_refseq[PROP] AND biomol_mrna[PROP]',
                   usehistory='y', retmax=0)
r = Entrez.read(h); h.close()
webenv, query_key, count = r['WebEnv'], r['QueryKey'], int(r['Count'])
print(f'{count} mRNAs queued on history server; use webenv/query_key with efetch')
```

### Inspect the translation before trusting a query

**Goal:** Catch translator misinterpretation before producing publication results.

**Approach:** Always print `QueryTranslation` for new queries and lock the rewritten string into the codebase as the canonical query.

```python
h = Entrez.esearch(db='pubmed', term='covid vaccine efficacy 2024', retmax=0)
r = Entrez.read(h); h.close()
print(r['QueryTranslation'])
# '("covid 19 vaccines"[MeSH Terms] OR ("covid 19"[All Fields] AND ...
# Now use this rewritten string explicitly to guarantee reproducibility.
```

### Discover fields for a database

```python
def list_fields(db):
    h = Entrez.einfo(db=db); r = Entrez.read(h); h.close()
    # Biopython 1.88 wraps DbInfo as a one-element list, not a dict -- index [0] first.
    return [(f['Name'], f['FullName'], f['Description']) for f in r['DbInfo'][0]['FieldList']]
```

### Spell-check before searching (catches typo-driven empty results)

```python
h = Entrez.espell(db='pubmed', term='breast canser')
r = Entrez.read(h); h.close()
print(r['CorrectedQuery'])  # 'breast cancer'
```

## Failure modes

### Silent retmax cap
Trigger and mechanism: see "retmax silent caps" above (9,999-record non-history cap).
- **Fix:** Always check `int(record['Count']) <= len(record['IdList'])`; switch to history server above ~5000.

### Query translation mismatch
- **Trigger:** Unqualified ambiguous term (e.g. `MARCH1` — Excel-renamed gene vs month abbreviation).
- **Mechanism:** EQT falls back to `[All Fields]` when no unambiguous mapping is found.
- **Symptom:** Either zero hits (gene symbol not in `[All Fields]`) or huge non-specific hits.
- **Fix:** Use field-qualified terms; for gene symbols, use HGNC ID via `gene` db lookup first.

### WebEnv expiration mid-pipeline
Trigger and mechanism: see "History server (WebEnv/QueryKey) semantics" above (TTL, idle eviction, error body).
- **Fix:** Parse error bodies (not just status codes); re-run ESearch and resume from `retstart`.

### Index lag for fresh deposits
Trigger and mechanism: see "Index lag" above.
- **Fix:** If the accession is known, use EFetch directly; only use ESearch for content-based discovery.

### Organism over-expansion
Trigger and mechanism: see "Organism field gotcha" above.
- **Fix:** Use `[Organism:exp]` to disable the walk, or constrain to a specific species/genus.

### Empty IdList with no error
- **Trigger:** Misspelled field name (`[gene]` works; `[gene_name]` returns nothing).
- **Mechanism:** Unknown field is silently coerced to `[All Fields]` — but combined with `AND` of a real field, the AND prunes everything.
- **Symptom:** Query that "should" match gets 0 results.
- **Fix:** Run EInfo on the db first to confirm field names; check `QueryTranslation`.

## Rate-limit math

| Auth | req/sec allowed | Sleep between calls | Bulk-friendly? |
|---|---|---|---|
| Email only | 3 | 0.34 s | Use history server, not parallel calls |
| Email + API key | 10 | 0.10 s | Modest parallelism (4 workers) is safe |
| Institutional bulk | Email `eutilities@ncbi.nlm.nih.gov` | Negotiated | For >100K queries; courtesy expected |

NCBI's terms of use ask that heavy automated queries run outside US weekday business hours (9 AM-5 PM ET). For genuinely bulk work, prefer the history server over parallel API calls — chunking against one session is faster and friendlier than scaling out.

## Common errors

| Error / symptom | Cause | Solution |
|---|---|---|
| `HTTPError 429` | Rate limit exceeded | Add `time.sleep(0.34)` or use API key |
| `HTTPError 400` | Field name or bracket malformed | Inspect EInfo field list; check brackets |
| `RuntimeError: ... email` | Missing `Entrez.email` | Set globally before any call |
| Empty `IdList`, large `Count` | Hit retmax cap | Set `retmax` explicitly or use history |
| `<ERROR>WebEnv not found</ERROR>` (HTTP 200) | Session expired | Re-run ESearch; parse XML body for errors |
| Query gives wildly wrong count | EQT misinterpretation | Print `QueryTranslation`; use field-qualified terms |

## References

- Sayers EW et al. (2024) Database resources of the National Center for Biotechnology Information in 2024. *Nucleic Acids Res* 52:D33-D43.
- Kans J. (2024) Entrez Direct: E-utilities on the Unix Command Line. NCBI Bookshelf NBK179288.
- NCBI. E-utilities In-Depth: Parameters, Syntax and More. NBK25499 (online manual; check current revision).
- Cock PJ et al. (2009) Biopython: freely available Python tools for computational molecular biology and bioinformatics. *Bioinformatics* 25:1422-1423.

## Related Skills

- entrez-fetch - Retrieve actual records once UIDs are in hand
- entrez-link - Cross-database navigation via ELink
- batch-downloads - History-server batch retrieval pipelines
- ncbi-datasets-cli - Modern alternative for genome/gene metadata (Datasets v2 CLI)
- geo-data - Specialized search semantics for the gds database
- sra-data - SRA metadata search before downloading runs
- ensembl-rest - Ensembl REST as alternative for Ensembl-native queries
