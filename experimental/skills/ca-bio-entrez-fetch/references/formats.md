# Formats



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

*Source-tested nucleotide ESummary docsums carry `AccessionVersion`, `Length`, `Title` and `TaxId`. Resolve `TaxId` via taxonomy for an organism name; interpreting the title's leading words is only a heuristic.

ESummary's documented hard limit is 10,000 docsums per call, but the practical sweet spot is ~500 (keeps the URL under length limits when IDs are comma-joined; for >500 use EPost to push IDs server-side first). Per-record payload is much smaller than EFetch. Use ESummary as the default for any metadata-only workflow.

## XML schema brittleness

`Entrez.read()` parses INSDSeq XML, PubmedArticle XML, Gene XML, etc. The schemas are NOT versioned; NCBI adds and renames fields without notice. Real-world consequence: a parser that worked in 2022 may KeyError in 2026 because a nested field moved.

Defensive patterns:
- Use `.get(key, default)` not `[key]` for every nested field
- For sequence content, prefer `SeqIO.read()` over `Entrez.read()` — the SeqIO parsers are versioned with BioPython
- Pin BioPython version in production code; expect to update the parser when NCBI changes the XML
- For PubMed, `Bio.Medline.parse(handle)` (against `rettype='medline'`) is more stable than the XML route



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
- **Fix:** Sniff the first line of the response with `expect_start()` (see references/setup.md) — `LOCUS` for GB, `>` for FASTA, `Run,` for SRA runinfo — and raise on mismatch. `Entrez.read()` routes already raise `RuntimeError` on a backend `<ERROR>` element.

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
