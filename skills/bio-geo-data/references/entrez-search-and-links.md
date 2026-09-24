# Entrez search and GEO links

Runnable versions live in `examples/`. Each sets `Entrez.email`, so edit it before running.

### Search GEO for studies matching a query

**Goal:** Find GSE accessions matching keywords + organism + study type.

**Approach:** ESearch on `gds` db with field-qualified terms; filter to `gse[Entry Type]`; summarize with ESummary.

`python examples/search_geo.py`: `search_geo(term, study_type, organism, gds_type, max_results)` builds the field-qualified query and returns the ESummary records; `detect_super_series(gse)` then reads the SOFT `!Series_relation` lines for each hit.

### Link GEO Series to SRA runs (preferred path: pysradb)

`python examples/geo_to_sra.py`: `gse_to_srr_pysradb(gse)` resolves GSE -> SRP -> SRR with `SRAweb().gse_to_srp()` then `srp_to_srr()`; `gse_to_srr_entrez(gse)` is the gds -> bioproject -> sra ELink fallback.

### Find datasets by PubMed citation

`python examples/geo_from_pubmed.py`: `find_geo_for_pubmed(pmid)` runs ELink `pubmed -> gds` then ESummary on the linked UIDs.
