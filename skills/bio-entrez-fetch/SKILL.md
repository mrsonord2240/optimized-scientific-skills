---
name: bio-entrez-fetch
description: Retrieve records from NCBI databases using Biopython Bio.Entrez (EFetch, ESummary). Use when downloading sequences, fetching GenBank/GenPept records, getting document summaries, parsing nested XML, navigating GI deprecation, choosing between rettype+retmode combinations, and parsing into Biopython SeqRecord/SwissProt objects. Covers nucleotide, protein, gene, pubmed, sra, gds, taxonomy, snp, clinvar.
tool_type: python
primary_tool: Bio.Entrez
license: MIT
author: GPTomics
---

## Version Compatibility

Reference examples tested with: BioPython 1.83+, Entrez Direct 21.0+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show biopython` then `help(Bio.Entrez.efetch)` to check signatures
- CLI: `efetch -version` then `efetch -help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Entrez Fetch

**"Download a record by accession from NCBI"** -> EFetch returns the full record content in a chosen format (FASTA, GenBank, XML, MEDLINE, etc.). ESummary returns a lightweight "docsum" object — much faster when only metadata is needed.

The agent's first decision is always: does this workflow need the full record, or just metadata? ESummary is 5-10x cheaper than EFetch for the equivalent record set. For "tell me the organism, length, and definition line for 10,000 accessions", ESummary wins by an order of magnitude.

This Skill assumes you already have UIDs or accessions to fetch. If you need to discover them from a search term, use `entrez-search` first (ESearch) — this Skill does not do discovery, except as the input side of the "History-server fetch" pattern (`references/history-server-fetch.md`).

- Python: `Entrez.efetch(db=..., id=..., rettype=..., retmode=...)` (BioPython)
- CLI: `efetch -db nucleotide -id NM_007294 -format gb` (Entrez Direct, NBK179288)
- R: `entrez_fetch(db=..., id=..., rettype=...)` (rentrez)

## Required Setup

```python
import time
from Bio import Entrez, SeqIO
Entrez.email = 'researcher@institution.edu'
# Entrez.api_key = '<your NCBI key>'  # optional: raises rate to 10 req/sec. A placeholder string here gives HTTP 400, so leave unset without a real key

def expect_start(text, prefix):
    '''EFetch can return an HTML/XML error body with HTTP 200 -- check the first line before parsing.'''
    if not text.lstrip().startswith(prefix):
        raise RuntimeError(f'Unexpected EFetch response, starts {text[:80]!r}')
    return text
```

`expect_start` is used by the plain-text patterns (`fetch_genbank` below, `sra_runinfo` in `references/sra.md`); see Failure modes -> "EFetch returns HTML error page".

## Decision matrix: rettype + retmode per database

The combinations are not orthogonal — each (db, rettype, retmode) triple is enabled or disabled by NCBI server-side. Wrong combinations return either silent empty responses or HTTP 400. The triples below are the safe, current set.

### nucleotide / protein

| rettype | retmode | Returns | Use when |
|---|---|---|---|
| `fasta` | `text` | FASTA | Just need sequence + defline |
| `gb` (nuc) / `gp` (prot) | `text` | Full flat file | Need annotations, features, references |
| `gbwithparts` | `text` | GB with CONTIG sequences inlined | Whole-genome shotgun assemblies; default `gb` returns CONTIG records requiring a chase to resolve |
| `fasta_cds_na` | `text` | CDS-only nucleotide | Extract coding regions from annotated GB |
| `fasta_cds_aa` | `text` | CDS-translated AA | Get translated proteins from GB record in one call |
| `xml` (== gb XML) | `xml` | INSDSeq XML | Programmatic parsing; the schema is unversioned and shifts |
| `acc` | `text` | Accession.version per line | Just resolve UID -> accession |
| `seqid` | `text` | Internal seq-id | Rarely needed |

### Other databases

| db | rettype / retmode | File |
|---|---|---|
| pubmed | `abstract`, `medline`, `xml` | `references/pubmed.md` |
| gene, taxonomy, gds | `gene_table`/`xml`, `xml`, default text | `references/gene-taxonomy-gds.md` |
| sra | `runinfo`, `xml` | `references/sra.md` |
| clinvar, snp | `vcv`/`xml`, `xml` | `references/variants-clinvar-snp.md` |

Wrong combinations return silent empty responses or HTTP 400; read the matching file before sending an unfamiliar one.

## GI deprecation (still bites in 2026)

NCBI stopped issuing new GI numbers for major nucleotide/protein submissions starting 2017. Records submitted after the cutoff have only `accession.version` identifiers. Many older scripts assume `id=<numeric_gi>`; passing a modern accession string also works, but mixing the two in one comma-separated id list is the bug.

**Rules:**
- For modern code, always pass `accession.version` strings.
- A bare accession without `.version` resolves to the latest version — fine for exploratory work, dangerous for reproducibility.
- Old `id=12345` GI lookups still work for records issued before 2017, but a search returning a UID that looks like a GI may actually be the legacy GI for an old record — assume UID is an opaque identifier.
- EFetch accepts comma-separated IDs of mixed types but the URL has a ~2000 char practical limit; chunk large ID lists into batches.

## ESummary vs EFetch triage

| Need | ESummary | EFetch (text) | EFetch (xml) |
|---|---|---|---|
| Title, length | yes | overkill | overkill |
| Organism | derived* | overkill | overkill |
| Authors of a PubMed article | yes | yes | yes |
| Full abstract text | no | `rettype=abstract` | better — structured |
| MeSH terms, grant info, PMC ID | no | no | yes |
| Sequence | no | `rettype=fasta` | overkill |
| Sequence features (CDS, exons) | no | `rettype=gb` | yes |
| Cross-references (xref) | partial | yes (in GB) | yes |
| Bulk metadata for 10K records | best (1 call per ~500) | slow | slow |

*Current nucleotide ESummary docsums (Biopython 1.88, checked 2026-09-17) carry no `Organism` field — only `AccessionVersion`, `Length`, `Title`, `TaxId`, etc. Derive organism from the leading binomial in `Title`, or resolve `TaxId` via a `db='taxonomy'` EFetch. See `bulk_summaries()` / `organism_of()` in `examples/fetch_summaries.py`.

ESummary's documented hard limit is 10,000 docsums per call, but the practical sweet spot is ~500 (keeps the URL under length limits when IDs are comma-joined; for >500 use EPost to push IDs server-side first). Per-record payload is much smaller than EFetch. Use ESummary as the default for any metadata-only workflow.

## XML schema brittleness

`Entrez.read()` parses INSDSeq XML, PubmedArticle XML, Gene XML, etc. The schemas are NOT versioned; NCBI adds and renames fields without notice. Real-world consequence: a parser that worked in 2022 may KeyError in 2026 because a nested field moved.

Defensive patterns:
- Use `.get(key, default)` not `[key]` for every nested field
- For sequence content, prefer `SeqIO.read()` over `Entrez.read()` — the SeqIO parsers are versioned with BioPython
- Pin BioPython version in production code; expect to update the parser when NCBI changes the XML
- For PubMed, `Bio.Medline.parse(handle)` (against `rettype='medline'`) is more stable than the XML route

## Code patterns

### Single sequence by accession

**Goal:** Fetch one nucleotide record as a SeqRecord with features.

**Approach:** EFetch with `rettype='gb', retmode='text'`; parse with `SeqIO.read()`.

**Reference (BioPython 1.83+):**
```python
from io import StringIO

def fetch_genbank(accession):
    h = Entrez.efetch(db='nucleotide', id=accession, rettype='gb', retmode='text')
    text = h.read(); h.close()
    return SeqIO.read(StringIO(expect_start(text, 'LOCUS')), 'genbank')

gb = fetch_genbank('NM_007294.4')
for feat in gb.features:
    if feat.type == 'CDS':
        print(feat.location, feat.qualifiers.get('product', ['?'])[0])
```

### Bulk metadata via ESummary

**Goal:** Get organism + length + title for 1,000 UIDs without downloading sequences.

**Approach:** ESummary on a comma-joined ID batch (max 500 per call by convention; supports 10K hard limit).

**Reference:** `examples/fetch_summaries.py` defines `bulk_summaries(db, ids, chunk=500)` (chunked ESummary, rate-limited) and `organism_of(docsum)` (derives the organism from `Title` when the docsum has no `Organism` field); run it for a demo, or copy the two functions.

### Extract CDS in one round-trip

**Goal:** Download the CDS-only translated protein sequences from a GenBank record without manually walking features.

**Approach:** Use `rettype='fasta_cds_aa'` — NCBI server-side extracts and translates every CDS in the record.

**Reference (BioPython 1.83+):**
```python
def cds_proteins(accession):
    h = Entrez.efetch(db='nucleotide', id=accession, rettype='fasta_cds_aa', retmode='text')
    return list(SeqIO.parse(h, 'fasta'))

proteins = cds_proteins('NC_000913.3')  # E. coli K-12 genome
print(f'{len(proteins)} CDS-translated proteins')
```

### Other code patterns

- PubMed with structured MeSH: `references/pubmed.md`
- History-server fetch (post-ESearch): `references/history-server-fetch.md`
- SRA UID -> SRR accession + run metrics: `references/sra.md`
- Taxonomy lineage by TXID: `references/gene-taxonomy-gds.md`
- ClinVar and dbSNP records: `references/variants-clinvar-snp.md`

## Failure modes

### Mixed-format batch silently truncates
- **Trigger:** Mixing modern accessions and legacy GIs in one comma-separated `id=`.
- **Mechanism:** EFetch parses left-to-right; on type-mismatch it may return only the prefix that succeeded.
- **Symptom:** Batch of 100 returns 47 records with no error.
- **Fix:** Validate that all IDs in a batch are the same type before sending.

### `gb` returns CONTIG instead of sequence
- **Trigger:** Fetching a whole-genome shotgun (WGS) assembly with `rettype='gb'`.
- **Mechanism:** Default GB output skips the contig sequence for assemblies, returning only the join() statement.
- **Symptom:** `len(record.seq) == 0` despite the record showing a length in metadata.
- **Fix:** Use `rettype='gbwithparts'` for assemblies; or for FASTA use `rettype='fasta'` directly.

### XML parse fails on schema drift
- **Trigger:** Code that was last touched in 2022 hits a new NCBI XML field layout.
- **Mechanism:** `Entrez.read()` uses cached DTDs that may not match current responses.
- **Symptom:** KeyError or ValidationError on a field that "always worked".
- **Fix:** Run `Entrez.read._XMLParser._DTDs.clear()` to force re-fetch of DTDs; upgrade BioPython; or switch to text format (`rettype='medline'` for pubmed) which is more stable.

### Silent empty response on bad rettype
- **Trigger:** Asking for `rettype='abstract'` on the nucleotide db (only valid for pubmed).
- **Mechanism:** EFetch returns empty text — no HTTP error.
- **Symptom:** `handle.read()` returns `''` or whitespace.
- **Fix:** Check the decision matrix above before sending unfamiliar combinations.

### Accession without `.version` returns wrong record later
- **Trigger:** Storing `'NM_007294'` (no version) for reproducibility years later.
- **Mechanism:** NCBI returns the current version, which may have changed annotation.
- **Symptom:** Re-run produces different CDS coordinates than the original analysis.
- **Fix:** Always pin `accession.version` (e.g. `NM_007294.4`); the version is in the GB LOCUS line.

### EFetch returns HTML error page
- **Trigger:** Invalid UID, mid-maintenance window, or expired WebEnv.
- **Mechanism:** Failure surfaces in HTML body, HTTP status is 200.
- **Symptom:** SeqIO chokes parsing HTML as GenBank; a naive CSV split (`sra_runinfo`) silently returns junk rows.
- **Fix:** Sniff the first line of the response with `expect_start()` (Required Setup) — `LOCUS` for GB, `>` for FASTA, `Run,` for SRA runinfo — and raise on mismatch. `Entrez.read()` routes already raise `RuntimeError` on a backend `<ERROR>` element.

## Common errors

| Error / symptom | Cause | Solution |
|---|---|---|
| `HTTPError 400` | Invalid id/db/rettype combo | Verify against decision matrix; check accession exists |
| `HTTPError 429` | Rate limit exceeded | Add `time.sleep(0.34)` or use API key |
| Empty `SeqRecord.seq` | WGS record with `rettype='gb'` | Use `rettype='gbwithparts'` |
| `ValueError: Sequence too short` | Wrong format declared to SeqIO | Match rettype to `SeqIO` format string |
| `ExpatError` | Got HTML where XML expected | Sniff response start; retry |
| KeyError on nested XML field | Schema drift | Use `.get()` defensively; pin BioPython |
| `TypeError: a bytes-like object is required, not 'str'` | Some databases (e.g. `sra`) return `bytes` from EFetch even with `retmode='text'` | `raw.decode() if isinstance(raw, bytes) else raw` before string ops |
| `TypeError: string indices must be integers` on a parsed XML element | Treating a `StringElement` (str subclass with `.attributes`) as a dict, e.g. `id['#text']` | Use `str(element)` for the value, `element.attributes.get(...)` for attributes |
| `ValueError: ... neither a DTD nor an XML Schema ...` | `Entrez.read()` on XML with no DTD (e.g. ClinVar VCV, some dbSNP responses) | Parse with `xml.etree.ElementTree` instead |


## Reference Files

| File | Read when |
|---|---|
| `references/pubmed.md` | Fetching PubMed abstracts, MEDLINE or XML (MeSH, grants, PMC ID) |
| `references/gene-taxonomy-gds.md` | Fetching `gene`, `taxonomy` lineage or `gds` records |
| `references/sra.md` | Converting SRA UIDs to SRR accessions and run metrics |
| `references/variants-clinvar-snp.md` | Fetching ClinVar VCV or dbSNP records |
| `references/history-server-fetch.md` | Pulling a large ESearch result set via `webenv`/`query_key` |

Runnable scripts (paths relative to this Skill folder): `scripts/history_fetch.py` (large ESearch result set to FASTA), `scripts/variant_records.py` (ClinVar / dbSNP record by UID). Runnable demos: `examples/`.

## References

- Sayers EW et al. (2024) Database resources of the National Center for Biotechnology Information in 2024. *Nucleic Acids Res* 52:D33-D43.
- Kans J. (2024) Entrez Direct: E-utilities on the Unix Command Line. NCBI Bookshelf NBK179288.
- NCBI. EFetch help. NBK25499.
- Cock PJ et al. (2009) Biopython: freely available Python tools for computational molecular biology and bioinformatics. *Bioinformatics* 25:1422-1423.

## Related Skills

- entrez-search - Find UIDs before fetching
- entrez-link - Cross-database navigation via ELink
- batch-downloads - History-server pipelines for large fetches
- ncbi-datasets-cli - Modern CLI for genome / gene metadata; often faster than EFetch
- sequence-io/read-sequences - Parse downloaded FASTA/GenBank with SeqIO
