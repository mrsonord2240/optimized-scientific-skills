# Formats



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
