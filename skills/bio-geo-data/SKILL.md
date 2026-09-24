---
name: bio-geo-data
description: Query and download from NCBI Gene Expression Omnibus (GEO) and EMBL-EBI's BioStudies/ArrayExpress mirror. Use when finding expression datasets, navigating SuperSeries vs SubSeries, choosing between series-matrix (submitter-normalized) and raw supplementary files, downloading via GEOparse (Python) or GEOquery (R/Bioconductor), linking GEO to SRA for raw reads, or distinguishing GSE/GSM/GPL/GDS record types. Encodes the SuperSeries trap, the series-matrix normalization-trust caveat, GEOmetadb deprecation, ArrayExpress migration to BioStudies, and processed-vs-raw decision matrix.
tool_type: mixed
primary_tool: Bio.Entrez
license: MIT
author: GPTomics
---

## Version Compatibility

Reference examples tested with: BioPython 1.83+, GEOparse 2.0+, R Bioconductor GEOquery 2.70+, pandas 2.2+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show biopython geoparse` then introspect signatures
- R: `packageVersion('GEOquery')`

Text encoding differs by language on Windows (checked 2026-09-21, Python 3.12 / R 4.4.3): Python's default follows the locale (`cp1252` here, unless `PYTHONUTF8=1`) while R's was UTF-8. GEO SOFT and matrix files contain UTF-8, so in Python always pass `encoding='utf-8', errors='replace'` to `gzip.open(..., 'rt')`, whatever the R side does.

If the GSE structure doesn't match expectations (missing fields, malformed series matrix), re-fetch from FTP directly and inspect the SOFT or MINiML file as source of truth.

# GEO Data

**"Pull expression data from GEO accession GSE..."** -> GEO stores Series (GSE), Samples (GSM), Platforms (GPL), and curated DataSets (GDS, frozen 2018). The single most consequential decision is **processed (series matrix) vs raw (supplementary files / linked SRA)** — the answer turns on how much trust the submitter's normalization deserves.

The single most-missed gotcha: **SuperSeries**. A GSE may be a meta-container (`!Series_relation = SuperSeries of: GSExxxxx`) holding multiple sub-studies on different platforms. Naively pulling samples from a SuperSeries gives mixed Affymetrix + Illumina + RNA-seq, mis-batched.

- Python: `Entrez.esearch(db='gds')`, GEOparse for full series download
- R: `GEOquery::getGEO()` (Bioconductor; more mature than GEOparse)
- CLI: `wget` from `ftp.ncbi.nlm.nih.gov/geo/series/...`

## Required Setup

```bash
pip install biopython GEOparse pandas pysradb
# OR for R-side:
# R: BiocManager::install('GEOquery')
```

```python
from Bio import Entrez
Entrez.email = 'researcher@institution.edu'
Entrez.api_key = 'optional'
```

## Workflow

1. Search gds db with field-qualified terms (`gse[Entry Type]`, `Homo sapiens[Organism]`, `expression profiling by high throughput sequencing[GDS Type]`); code in `references/entrez-search-and-links.md`.
2. For any GSE returned, check `!Series_relation` in SOFT to detect SuperSeries before pulling.
3. Pick the right download path: series matrix for fast-and-trusting; supplementary files for raw Affymetrix / submitter counts; SRA-link for RNA-seq raw FASTQ.
4. Read `!Sample_data_processing` to surface what's actually in the series matrix.
5. For R-side analyses, recommend GEOquery (Bioconductor) over GEOparse for supplementary file reliability; see `references/geoparse-geoquery.md`.
6. For SRA hand-off, use pysradb to resolve GSE -> SRP -> SRR (`references/entrez-search-and-links.md`); pass run list to sra-data skill.
7. Warn on stale GEOmetadb usage; recommend pysradb / Entrez gds (`references/legacy-and-formats.md`).
8. For ArrayExpress accessions (E-MTAB-*), use the new BioStudies URL (`references/legacy-and-formats.md`).

## GEO record taxonomy

| Prefix | Type | Granularity | What's in it |
|---|---|---|---|
| GSE | Series | One study | Title, summary, design, links to GSMs, supplementary files |
| GSM | Sample | One biological/technical sample | Submitter metadata, per-sample processed data, link to raw SRA |
| GPL | Platform | One array / sequencer | Probe annotations or sequencer model |
| GDS | DataSet | Curated, normalized subset of one GSE | Re-normalized expression matrix (frozen 2018; new GDS no longer created) |
| GSEXXX SuperSeries | Series meta-container | Wraps multiple SubSeries | `!Series_relation = SuperSeries of: ...` |

**`GDS` is dead-as-format**: NCBI stopped creating new GDS records in 2018. Existing GDS still queryable but use GSE for anything current.

**`GDS` is overloaded**: Entrez's `db='gds'` is the query endpoint indexing all four record types (GSE/GSM/GPL/GDS) despite the name; don't confuse it with the `GDS` record-type prefix above, which means specifically the frozen curated DataSet type.

## The SuperSeries trap

A SuperSeries (GSE) wraps multiple SubSeries, often with different platforms. Detection:

```python
# Read the !Series_relation field from SOFT format
from Bio import Entrez
h = Entrez.esummary(db='gds', id='200122288')   # example
r = Entrez.read(h)[0]; h.close()
print(r.get('summary'))   # may or may not flag SuperSeries
# Definitive check: download SOFT and grep:
#   curl ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE122nnn/GSE122288/soft/GSE122288_family.soft.gz | zgrep Series_relation
```

A `SuperSeries of: GSE12345` line means the SuperSeries' samples are the union of all SubSeries — almost certainly mixed-platform / mixed-batch. Process each SubSeries independently.

Symmetric trap: a paper may cite a SubSeries (`SubSeries of: GSEsuper`) where the wider context is essential — check both directions.

## Decision matrix: processed vs raw vs SRA

| Question | Source | Trust level |
|---|---|---|
| "I want expression values; submitter normalization is fine" | Series matrix (`GSE_series_matrix.txt.gz`) | Trust submitter's normalization |
| "I want raw Affymetrix CEL files and to do my own RMA" | Supplementary files (`suppl/`) | Re-normalize locally |
| "I want raw RNA-seq FASTQ" | pysradb `gse_to_srp -> srp_to_srr` (Entrez gds->sra ELink unreliable; code in `references/entrez-search-and-links.md`) | Always raw; processed at submitter is rarely re-usable |
| "I want submitter-provided counts (RNA-seq)" | Supplementary files (usually a `*_counts.txt.gz`) | Trust at risk; submitter pipelines vary |
| "I want a curated subset across many studies" | Use ArchS4 (https://archs4.org) or recount3 | Curated re-processing |

**Default to raw whenever possible.** For Affymetrix: CEL + locally-run RMA is far more reliable than the submitter's "normalized" matrix. For RNA-seq: SRA FASTQ + locally-run alignment/quantification is the only reproducible path; submitter counts often use a private pipeline.

**Determining platform technology (array vs RNA-seq):** an ESummary record's `gdsType` field states this directly (e.g. `"Expression profiling by high throughput sequencing"` -> RNA-seq; `"Expression profiling by array"` / `"Methylation profiling by array"` -> array-based, likely Affymetrix or Illumina):

```python
h = Entrez.esearch(db='gds', term='GSE147507[Accession]', retmax=1)
s = Entrez.read(h); h.close()
h = Entrez.esummary(db='gds', id=s['IdList'][0])
gse_record = Entrez.read(h)[0]; h.close()
is_sequencing = 'high throughput sequencing' in gse_record.get('gdsType', '')
```

## Series matrix files

A series matrix (`GSE12345_series_matrix.txt.gz`) is a header (sample metadata as `!Sample_*` lines) plus a sample-by-feature expression table. The format is fragile and the values' provenance is whatever the submitter chose. Critical caveats:

- For Affymetrix: the matrix is usually RMA-normalized but submitters sometimes apply additional transforms (log2, scaling, batch correction).
- For RNA-seq: the matrix is sometimes log-CPM, sometimes raw counts, sometimes VST/rlog — read `!Series_overall_design` and `!Sample_data_processing` to know.
- The header has `!Sample_characteristics_ch1` rows that hold the metadata of interest — these are submitter-formatted strings, often inconsistent within one series.
- `!Sample_data_processing` can be **entirely absent**, not just terse (e.g. GSE470 has zero such lines) — always read it via `metadata.get('!Sample_data_processing', [])`, never direct key access, and treat an empty result as "field not provided," not "no processing was done."

## Code patterns

### Detect SuperSeries before pulling data

**Goal:** Avoid mixing platforms by detecting SuperSeries structure first.

**Approach:** Stream the SOFT family file and read the `!Series_relation` keys in its header (stops at the first `^PLATFORM`; family files can be hundreds of MB).

```bash
python scripts/geo_series.py relation GSE346738
# {'super_of': ['GSE283260', 'GSE346737'], 'sub_of': None}  -> SuperSeries; process subseries separately
# (checked live 2026-09-21; SuperSeries status drifts as submitters restructure series, re-verify before relying on a fixed example)
```

As a library: `from geo_series import check_super_or_sub_series` (returns `{'super_of': [...], 'sub_of': ...}`).

### Download series matrix with submitter caveat

```bash
python scripts/geo_series.py matrix GSE470
# GSE470: 12625 features x 12 samples; then every distinct !Sample_data_processing note (none for GSE470)
python scripts/geo_series.py selftest   # offline regression test: non-ASCII fixtures, encoding= on every gzip.open 'rt'
```

As a library: `from geo_series import download_series_matrix, parse_series_matrix`; `parse_series_matrix(path)` returns `(metadata, expr)`, where `metadata` maps each `!` header key to its list of values.

## Reference Files

| File | Read when |
|---|---|
| `references/entrez-search-and-links.md` | Searching `gds` with field-qualified terms, resolving GSE -> SRP -> SRR with pysradb, or finding GEO datasets from a PMID |
| `references/geoparse-geoquery.md` | Choosing between GEOparse and GEOquery, or downloading a full Series or supplementary files in Python or R |
| `references/legacy-and-formats.md` | Choosing SOFT vs MINiML, handling GEOmetadb pipelines, or old ArrayExpress (E-MTAB-*) URLs |

## Failure modes

### SuperSeries pulled as one experiment
- **Trigger:** GSE accession from a paper; turns out to be a SuperSeries wrapping multiple platforms.
- **Mechanism:** Default download merges all samples without flagging the structure.
- **Symptom:** Downstream batch correction can't recover the mixed-platform structure; spurious "batch" effects.
- **Fix:** Always check `!Series_relation` in SOFT before pulling; process SubSeries independently.

### Series matrix is not what it appears to be
- **Trigger:** Series matrix downloaded; treated as RMA-normalized when submitter applied additional transforms.
- **Mechanism:** Series matrix contents are at submitter's discretion.
- **Symptom:** Re-analysis gives different answers than the published paper.
- **Fix:** Read `!Sample_data_processing` to know what's in the matrix; re-normalize from raw if in doubt.

### Submitter-provided RNA-seq counts mis-trusted
- **Trigger:** Using a `*_counts.txt.gz` supplementary file as the count matrix.
- **Mechanism:** Submitter's pipeline (aligner, GTF version, counting strategy) is rarely documented.
- **Symptom:** Counts don't agree with re-quantification from SRA FASTQ.
- **Fix:** Pull SRA FASTQ + re-quantify with a known pipeline (Salmon, kallisto, STAR + featureCounts).

### Platform GPL mismatch
- **Trigger:** One GSE with multiple platforms; series matrix split across multiple files.
- **Mechanism:** `GSE_series_matrix.txt.gz` is the merged one; per-platform are `GSE-GPLxxx_series_matrix.txt.gz`.
- **Symptom:** "Missing samples" or NaN-heavy expression matrix.
- **Fix:** Download per-platform matrix files; check `!Series_platform_id` count.

### GEOparse supplementary files flakey
- **Trigger:** `gse.download_supplementary_files()` silently misses files.
- **Mechanism:** Known issue with the GEOparse FTP enumeration since ~2022.
- **Symptom:** Local cache missing CEL or counts files.
- **Fix:** Use R GEOquery or direct FTP `wget -r` on the suppl/ subdirectory.

### ArrayExpress URL rot
- **Trigger:** Old paper links `https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-1234/`.
- **Mechanism:** ArrayExpress migrated to BioStudies in 2020.
- **Symptom:** 404 or redirect.
- **Fix:** Use `https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-1234`.

### GEOmetadb stale
- **Trigger:** Old pipeline downloads `GEOmetadb.sqlite` for fast queries.
- **Mechanism:** GEOmetadb unmaintained since 2020.
- **Symptom:** Missing recent series; outdated annotations.
- **Fix:** Switch to pysradb for SRA-linked queries; Entrez gds for full GEO.

## Common errors

| Error / symptom | Cause | Solution |
|---|---|---|
| Empty IdList for `gse[entry_type]` | Wrong field name | Use `gse[Entry Type]` (case-sensitive) |
| Matrix file has no expression data | SuperSeries with no aggregate matrix | Pull per-SubSeries matrices |

## References

- Edgar R, Domrachev M, Lash AE. (2002) Gene Expression Omnibus: NCBI gene expression and hybridization array data repository. *Nucleic Acids Res* 30:207-210.
- Barrett T, Wilhite SE, Ledoux P, et al. (2013) NCBI GEO: archive for functional genomics data sets - update. *Nucleic Acids Res* 41:D991-D995.
- Davis S, Meltzer PS. (2007) GEOquery: a bridge between the Gene Expression Omnibus (GEO) and BioConductor. *Bioinformatics* 23:1846-1847.
- Gumienny R. GEOparse: Python library to parse GEO databases. https://github.com/guma44/GEOparse (no journal publication).
- Sarkans U, Gostev M, Athar A, et al. (2018) The BioStudies database--one stop shop for all data supporting a life sciences study. *Nucleic Acids Res* 46:D1266-D1270.
- Lachmann A, Torre D, Keenan AB, et al. (2018) Massive mining of publicly available RNA-seq data from human and mouse. *Nat Commun* 9:1366. (ARCHS4)
- Wilks C, Zheng SC, Chen FY, et al. (2021) recount3: summaries and queries for large-scale RNA-seq expression and splicing. *Genome Biol* 22:323.

## Related Skills

- entrez-search - General gds search
- entrez-link - gds <-> pubmed, bioproject links (gds->sra ELink is unreliable; use pysradb)
- sra-data - Download raw FASTQ from GEO-linked SRA runs
- expression-matrix/normalization - Re-normalize raw expression data
- rna-quantification/alignment-free-quant - Salmon/kallisto re-quantification of GEO/SRA data
- ensembl-rest - Cross-reference Ensembl IDs in series-matrix files
