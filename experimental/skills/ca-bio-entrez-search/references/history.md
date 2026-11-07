# History



| Endpoint behavior | Cap | Workaround |
|---|---|---|
| Default `retmax` | 20 | Set explicitly |
| Legacy esearch.fcgi (no `usehistory`) | **9,999** silent cap | Use history server |
| `usehistory='y'` + ESearch | 100,000 per page | Page with `retstart` against the WebEnv |
| EPost (to push IDs server-side) | No small universal cap | Batch a few hundred IDs per call and union with QueryKey |

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

Entrez and PubMed indexes update on provider-controlled schedules, so a fresh deposit or newly assigned
MeSH term can be absent temporarily. Do not infer that a known accession is invalid from an empty
search result: use EFetch for the known accession, or NCBI Datasets for genome metadata, and record
the EInfo `LastUpdate` value when timing matters.



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
