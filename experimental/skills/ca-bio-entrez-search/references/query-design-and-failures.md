# Query Design and Failure Recovery

Read this reference when selecting database-specific fields, validating an unexpected count, or
recovering from an Entrez request failure. It supplements the operational workflow in `SKILL.md`.

## Field-qualified query patterns

| Database | Common fields | Notes |
|---|---|---|
| pubmed | `[Title]`, `[TIAB]`, `[MeSH]`, `[Author]`, `[Journal]`, `[PDAT]`, `[DCOM]`, `[PMC]` | `[TIAB]` is more permissive than `[Title]`; `[MeSH]` can lag. `pubmed pmc[sb]` selects the PMC full-text subset, not a separate database. |
| nucleotide | `[Organism]`, `[Gene Name]`, `[Accn]`, `[SLEN]`, `[Filter]`, `[PROP]` | `srcdb_refseq[PROP]` limits to RefSeq; `biomol_genomic[PROP]` filters molecule type. |
| protein | `[Organism]`, `[Gene Name]`, `[Accn]`, `[MOLWT]`, `[PROP]` | `swissprot[Filter]` limits to reviewed records. |
| gene | `[Gene/Locus]`, `[Organism]`, `[Chromosome]`, `[Gene Type]` | `[Gene Type]` includes `protein-coding`, `pseudo`, and `ncRNA`. |
| sra | `[Organism]`, `[Platform]`, `[Strategy]`, `[Library Source]`, `[BioProject]` | `[Strategy]` accepts `RNA-Seq`, `WGS`, `ChIP-Seq`, and related values. |
| gds (GEO) | `[Organism]`, `[Entry Type]`, `[GDS Type]`, `[Platform]` | `gse[Entry Type]` selects Series; `gds[Entry Type]` selects curated DataSets. |
| taxonomy | `[Scientific Name]`, `[Common Name]`, `[Rank]`, `[TXID]` | TXID is the numeric taxonomy ID. |
| clinvar | `[Gene Name]`, `[Clinical Significance]`, `[Variation Type]` | Use `pathogenic[CLIN]` for pathogenic records only. |

```python
# Curated human RefSeq mRNAs, 500--5000 nt
term = 'Homo sapiens[ORGN] AND srcdb_refseq[PROP] AND biomol_mrna[PROP] AND 500:5000[SLEN]'

# Reviewed human kinases
term = 'Homo sapiens[ORGN] AND swissprot[Filter] AND kinase[Protein Name]'

# Recent human studies with full text in PMC
term = 'CRISPR[Title] AND humans[MeSH Terms] AND last 30 days[EDAT] AND pubmed pmc[sb]'
```

`[Organism]` (and `[ORGN]`) is taxonomy-walked: `mammalia[ORGN]` includes descendant species.
Use `[Organism:noexp]` for records tagged at exactly that node; `:exp` explicitly retains the
default expansion. First verify the exact field names
with EInfo, then inspect `QueryTranslation`; unknown or ambiguous fields may be coerced to
`[All Fields]` rather than failing loudly.

## Failure modes and recovery

| Symptom | Likely cause | Recovery |
|---|---|---|
| `Count` exceeds returned `IdList` length | `retmax` truncation; non-history searches can silently cap at 9,999 | Set `retmax` explicitly, then switch to the history server for large sets. |
| Huge, zero, or surprising count | Entrez Query Translator reinterpreted an ambiguous or unqualified term | Print `QueryTranslation`, use field-qualified terms, and resolve ambiguous gene symbols through `gene`. |
| `<ERROR>WebEnv not found</ERROR>` with HTTP 200 | History session expired | Parse the body, re-run ESearch, and resume at the saved `retstart`. |
| Freshly deposited accession is absent | Search index or MeSH indexing lag | Use EFetch with the known accession; use ESearch for content discovery. |
| Taxonomic query expands unexpectedly | `[Organism]` includes descendants | Use `[Organism:exp]` or a narrower taxon. |
| `HTTPError 429` | Rate limit exceeded | Sleep 0.34 s without an API key, 0.10 s with one; do not parallelize one history session. |
| `HTTPError 400` | Malformed field/brackets or excessively long term | Check EInfo field names and brackets; validate the term, then EPost IDs in batches. |
| Biopython warning that email is unspecified | `Entrez.email` was not set | Set it globally before the first Entrez call; NCBI uses it to contact a high-volume caller. |
