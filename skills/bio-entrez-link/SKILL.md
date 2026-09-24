---
name: bio-entrez-link
description: Find cross-database references between NCBI databases using Biopython Bio.Entrez (ELink). Use when navigating gene to protein/structure, sequence to publication, PubMed to GEO, BioProject to SRA runs, or discovering all link relationships for a record. Covers linkname semantics, cmd= variants, asymmetric link warnings, neighbor_history for >200 input IDs, and per-database link tables.
tool_type: python
primary_tool: Bio.Entrez
license: MIT
author: GPTomics
---

## Version Compatibility

Reference examples tested with: BioPython 1.83+, Entrez Direct 21.0+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show biopython` then `help(Bio.Entrez.elink)` to check signatures
- CLI: `elink -version` then `elink -help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Entrez Link

**"Find records linked to this record in another NCBI database"** -> ELink walks the curated, weekly-maintained link tables between Entrez databases. A link is an asserted relationship (e.g. "this PubMed article describes this nucleotide sequence"), not a similarity hit.

ELink is the navigation layer of Entrez. The decision that matters most is **which `linkname` to use** — not which databases. A single (`dbfrom`, `db`) pair can have a dozen `linkname` variants distinguishing curation level, evidence type, and direction. Picking the wrong one changes the result set by a multiple (BRCA1, live 2026-09-21: 368 `gene_protein_refseq` proteins vs 1087 `gene_protein`).

- Python: `Entrez.elink(dbfrom=..., db=..., id=..., linkname=...)` (BioPython)
- CLI: `elink -db pubmed -target gene -name pubmed_gene_rif` (Entrez Direct)
- R: `entrez_link(dbfrom=..., db=..., id=...)` (rentrez)

## Required Setup

```bash
pip install biopython
```

```python
from Bio import Entrez
Entrez.email = 'researcher@institution.edu'
Entrez.api_key = 'optional_api_key'  # raises rate to 10 req/sec
```

## Workflow

1. Set `Entrez.email` and `Entrez.api_key`.
2. Confirm each source UID resolves in `dbfrom` and is the record you meant: `python scripts/check_source_ids.py <dbfrom> <uid> ...` (see Failure modes, "Mismatched dbfrom and id namespace").
3. For unfamiliar (`dbfrom`, `db`) pairs, run `cmd='acheck'` first to enumerate linknames.
4. Pick a `linkname` deliberately — prefer curated variants (`*_refseq`, `*_rif`, `*_swissprot`) for analyses; use the umbrella `gene_protein` only for exploration. Per-pair linknames for other databases: `references/link_catalog.md`.
5. For >200 source IDs, EPost first and ELink with `cmd='neighbor_history'` (code: `references/code_patterns.md`, or `examples/chain_links.py`).
6. Pass source IDs as a Python **list** to get one LinkSet per input UID (a comma-joined string returns a single merged LinkSet, the union), and iterate the response — never assume a single LinkSetDb covers all inputs.
7. Guard for empty `LinkSetDb` before indexing.
8. Document the asymmetry of round-trip queries when results matter for publication.

## The `linkname` decision (most important)

For most (`dbfrom`, `db`) pairs NCBI exposes multiple link tables. The qualifiers in the name encode the curation level and the evidence source. Choose deliberately.

### gene -> protein (representative example)

| linkname | Returns | When to use |
|---|---|---|
| `gene_protein` | All linked proteins (curated + automated) | Exploration; expect several times more hits (BRCA1: 1087 vs 368 RefSeq, TP53: RefSeq 25) |
| `gene_protein_refseq` | RefSeq proteins only (every RefSeq isoform, not one canonical protein) | Reference-quality analyses; orthology |
| `gene_protein_swissprot` | Reviewed UniProt entries with NCBI cross-ref | Functional annotation; literature support |

### pubmed -> gene

| linkname | Returns |
|---|---|
| `pubmed_gene` | Genes mentioned in this paper (text-mined + curated) |
| `pubmed_gene_rif` | Genes with a Reference Into Function (curated, high-quality) |
| `pubmed_gene_pubmed` | Other PubMed records sharing gene linkage (rare use) |

### nucleotide -> protein

| linkname | Returns |
|---|---|
| `nuccore_protein` | All proteins encoded by this nucleotide record (CDS-linked) |
| `nuccore_protein_refseq` | RefSeq proteins only |

### Discover what link names exist for a pair

```python
h = Entrez.elink(dbfrom='gene', db='protein', id='672', cmd='acheck')
record = Entrez.read(h); h.close()
for ls in record[0]['IdCheckList']['IdLinkSet'][0]['LinkInfo']:
    print(f'{ls["LinkName"]}  -> {ls["DbTo"]} | {ls.get("MenuTag", "<none>")} ({ls.get("HtmlTag", "<none>")})')
```

`cmd='acheck'` is the only authoritative way to enumerate available linknames — they change with each NCBI release. LinkInfo entries are keyed `LinkName` (not `Name`) on current Biopython/NCBI, and some entries omit `MenuTag`/`HtmlTag` entirely (e.g. `nuccore_nuccore_mrnaonly`, `pubmed_pmc_local`) — always access both with `.get(..., '<none>')`, never index directly.

## Decision table: which `cmd` for which goal

| Goal | cmd | Returns |
|---|---|---|
| Get linked records | `neighbor` (default) | Linked IDs in target db |
| Get linked + relevance scores | `neighbor_score` | IDs with similarity scores (mostly `pubmed_pubmed`) |
| Get >200 source IDs in one go | `neighbor_history` | WebEnv + QueryKey for downstream EFetch |
| Enumerate available links | `acheck` | List of all linknames for source IDs |
| Check if any link exists | `ncheck` | Boolean per source ID |
| Check specific link exists | `lcheck` | Boolean per source ID + linkname |
| Get NCBI HTML link URLs | `llinks` | URLs to Entrez record pages |
| Get external provider links | `prlinks` | URLs to journal sites, etc. |

The `neighbor_history` cmd is essential when source `id` count exceeds ~200 — past that, the URL-length limit makes the comma-joined form fail. With `neighbor_history` ELink puts results on the history server and returns WebEnv/QueryKey for downstream pickup.

## Asymmetric link warning

ELink relationships are **not guaranteed symmetric**. `pubmed_gene` and `gene_pubmed` may return different sets because:
- Direction-dependent curation: gene-to-PubMed is curated by NCBI staff (GeneRIF); PubMed-to-gene includes text-mining.
- Cutoffs: some link tables truncate at N best links in one direction but not the other.
- Index lag asymmetry: when one db updates faster than the other.

If round-trip consistency matters (e.g. "every gene mentioned in this paper, then every paper mentioning each gene"), expect the round-trip set to be larger than the input — and never assume `A -> B -> A` returns the original ID alone. The original PMID may be missing from the round-trip set while new PMIDs appear. Document the directional asymmetry, and use the more-curated linkname (`*_rif` variants) when fidelity matters.

## Clinically-actionable link tables

`gene_clinvar`, `gene_omim`, `gene_gtr`, and `gene_medgen_diseases` return real, curated clinical-variant and disease data — not a toy example (e.g. `gene_clinvar` on BRCA1 alone returns 16,000+ linked ClinVar records). When one of these link types is surfaced in response to a request framed around a specific patient or personal diagnosis, do not make a diagnostic or prescriptive claim from the link count or record set — include a clinician/genetic-counselor referral instead. This applies regardless of general safety training: the caveat is part of correct output for these linknames specifically.

## Reference Files

| File | Read when |
|---|---|
| `references/link_catalog.md` | You need the linknames for a specific (`dbfrom`, `db`) pair (gene, nuccore, protein, pubmed, bioproject), beyond the gene/pubmed/nucleotide tables above or `cmd='acheck'` |
| `references/code_patterns.md` | You are writing ELink code: single/batch/history-server linking, chaining, `neighbor_score`, BioProject -> SRA |

## Failure modes

### Wrong linkname multiplies the result set
- **Trigger:** Using `gene_protein` when `gene_protein_refseq` was intended.
- **Mechanism:** `gene_protein` includes all automated and predicted entries (XP_* RefSeq plus all GenBank submissions).
- **Symptom:** Several times more proteins than intended (BRCA1: 1087 instead of 368 RefSeq isoforms), including GenBank-submission and predicted entries.
- **Fix:** Pick the curated linkname; verify counts on a known gene. `_refseq` still returns every RefSeq isoform (368 for BRCA1), so filter to one isoform (e.g. MANE Select) if a single canonical protein is wanted.

### Empty LinkSetDb on valid input
- **Trigger:** Gene with no linked records in the requested target.
- **Mechanism:** `record[0]['LinkSetDb']` is an empty list, not raising an error.
- **Symptom:** `KeyError` if code assumes `record[0]['LinkSetDb'][0]` always exists.
- **Fix:** Always guard `if not record[0]['LinkSetDb']: return []`.

### URL length limit on large batches
- **Trigger:** Comma-joined `id=` with 200+ IDs.
- **Mechanism:** HTTP GET URL exceeds NCBI's parsing limit (~2000 chars).
- **Symptom:** HTTP 414 URI Too Long, or silent truncation.
- **Fix:** EPost the IDs first, then ELink with `cmd='neighbor_history'`.

### Chunked EPost links only the last chunk
- **Trigger:** EPosting >200 IDs in several calls into one WebEnv, then ELinking with the last `QueryKey`.
- **Mechanism:** Each EPost creates its own QueryKey holding only that chunk; posting into an existing WebEnv does not merge. Live 2026-09-21: 250 gene UIDs posted as 200 + 50, ELink from the last key returned 197 proteins, from the union 1374.
- **Symptom:** A valid WebEnv/QueryKey and no error, but the linked set covers only the final chunk.
- **Fix:** Union the chunk keys before linking: `Entrez.esearch(db=dbfrom, term='#1 OR #2', WebEnv=webenv, usehistory='y', retmax=0)` and link from the returned `QueryKey` (implemented in `examples/chain_links.py`, `link_batch_via_history`). Check the set size with `esearch(db=target, term='#<key>', WebEnv=..., retmax=0)['Count']`.

### Comma-joined vs list `id`, indexing confusion
- **Trigger:** Sending several IDs, then reading `record[0]['LinkSetDb'][0]['Link']` as if it covered all of them, or mapping results back to inputs after a comma-joined call.
- **Mechanism:** `id=['672','7157']` (list) sends one `id=` per UID and returns one `LinkSet` per input, in input order. `id='672,7157'` (comma-joined string) returns a single `LinkSet` whose `IdList` holds both UIDs and whose links are the union. Live 2026-09-21: list -> 2 linksets (368 and 25 proteins); comma-joined -> 1 linkset (393).
- **Symptom:** List form read as one set: only the first input's links are processed, the rest are dropped. Comma-joined form read per input: `linkset['IdList'][0]` is only the first UID and the others' links are attributed to it.
- **Fix:** For per-input results pass a list, iterate `for linkset in record:` and map by `linkset['IdList'][0]`. Use the comma-joined form only when the union is what you want.

### Invalid linkname (HTTP 400)
- **Trigger:** A `linkname` that does not exist for the (`dbfrom`, `db`) pair.
- **Symptom:** `HTTPError 400`.
- **Fix:** Enumerate valid linknames with `cmd='acheck'`.

### Mismatched dbfrom and id namespace
- **Trigger:** Passing a PMID into `dbfrom='nucleotide'`, or a nucleotide UID into `dbfrom='pubmed'` (a common slip when chaining UIDs across steps).
- **Mechanism:** ELink returns no error (`ERROR: []`) — it looks the number up in `dbfrom`, finds nothing, and returns the same empty `LinkSetDb` as a genuine "no links" answer. Numeric UIDs also collide across databases (PubMed 31322957 and nucleotide 31322957 are unrelated records), so a wrong `dbfrom` can silently link the wrong record.
- **Symptom:** Empty LinkSetDb on a "valid" ID, or links for a record you did not mean.
- **Fix:** Before linking, run `python scripts/check_source_ids.py <dbfrom> <uid> ...` (or `checked_elink()` from the same file): it ESummary-resolves every UID in `dbfrom`, prints its title/caption for you to compare with the record you meant, and exits non-zero on a UID that does not resolve (PMIDs are db=pubmed, GeneIDs are db=gene; accessions and symbols are not UIDs — resolve them with entrez-search first).

## References

- Sayers EW et al. (2024) Database resources of the National Center for Biotechnology Information in 2024. *Nucleic Acids Res* 52:D33-D43.
- Kans J. (2024) Entrez Direct: E-utilities on the Unix Command Line. NCBI Bookshelf NBK179288.
- NCBI. ELink help. NBK25499.

## Related Skills

- entrez-search - Resolve UIDs before linking
- entrez-fetch - Retrieve linked records' content
- batch-downloads - History-server retrieval after ELink with `neighbor_history`
- geo-data - Specialized gds <-> pubmed/bioproject links (gds->sra ELink unreliable; use pysradb)
- ncbi-datasets-cli - Modern alternative for gene/genome cross-reference queries
- local-blast / remote-homology - ELink returns asserted database relationships, not a similarity hit; for actual sequence/structure similarity, use these instead
