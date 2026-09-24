---
name: bio-sra-data
description: Download raw sequencing reads from NCBI SRA using sra-tools (prefetch, fasterq-dump, vdb-validate) or the ENA mirror. Use when pulling FASTQ for SRR/ERR/DRR accessions, deciding between SRA-direct, ENA mirror, or AWS/GCP cloud mirror (STRIDES), handling --include-technical for 10x and other single-cell records, validating with MD5/vdb-validate, navigating SRR/SRX/SRS/SRP/PRJNA hierarchy, or finding accessions via pysradb. Encodes SRA cloud-egress economics, the fasterq-dump uncompressed-scratch trap, and the --max-size default that silently truncates large prefetches.
tool_type: cli
primary_tool: sra-tools
license: MIT
author: GPTomics
---

## Version Compatibility

Reference examples tested with: sra-tools 3.0+ (fasterq-dump, prefetch, vdb-validate, vdb-config), pysradb 2.2+, ENA portal API 2.0+

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `fasterq-dump --version`, `prefetch --version`
- Python: `pip show pysradb`

If a flag is unrecognized or behavior changes, run `<tool> --help` and adapt.

# SRA Data

**"Download FASTQ from this SRA accession"** -> Two paths exist in 2026: the **SRA toolkit** (NCBI's official, with prefetch + fasterq-dump) and the **ENA mirror** (EMBL-EBI's mirror with direct FASTQ download, often faster). For >1 TB workflows, a third path: **AWS Open Data** (STRIDES program) where same-region EC2 pulls SRA data with zero egress cost.

The single most impactful decision is **where to pull from**. SRA-direct is the default but ENA is faster more often than not, and AWS Open Data is the right answer for cloud-native analysis pipelines.

- CLI: `prefetch SRR...`, `fasterq-dump SRR...`, `vdb-validate SRR...` (sra-tools)
- CLI: `curl https://ftp.sra.ebi.ac.uk/...` (ENA mirror; direct FASTQ)
- CLI: `aws s3 cp s3://sra-pub-run-odp/sra/SRR.../SRR... ./SRR....sra ...` (STRIDES; object is unsuffixed; same-region free)
- Python: `pysradb` for metadata; `subprocess` for download

## Required Setup

```bash
# sra-tools (toolkit)
conda install -c bioconda sra-tools           # 3.0+
fasterq-dump --version                        # confirm

# Configure cache location (default ~/ncbi/ -- often too small)
vdb-config --cfg                              # show current config
vdb-config --set /repository/user/main/public/root=/data/sra_cache

# Optional: pysradb for metadata
pip install pysradb
```

For STRIDES cloud, install `aws-cli` (no NCBI auth needed for public buckets) and run from EC2 in the bucket's own region (`us-east-1` for AWS, `us-central1` for GCP) for free egress; a cross-region pull is charged. Commands are under "Cloud (STRIDES) access".

## Decision matrix: where to pull from

| Source | When best | Speed | Cost |
|---|---|---|---|
| **ENA mirror** (FTP/Aspera) | Default for most workflows | Often fastest; direct FASTQ (no SRA->FASTQ conversion needed) | Free; no rate limit observed |
| **SRA toolkit + AWS STRIDES** | Same-region EC2/EKS | Fastest within AWS us-east-1 | Free egress within region; small storage cost |
| **SRA toolkit + GCP STRIDES** | Same-region GCP Compute Engine | Fastest within GCP us-central1 | Free egress within region |
| **SRA-direct (prefetch + fasterq-dump)** | On-prem; small downloads; need SRA-format access | Variable; can be slow off-peak fails | Free; NCBI throttles by IP |
| **Aspera (`ascp`)** | Institutional accounts only | Faster than HTTPS on long links | NCBI public Aspera retired 2019; ENA public Aspera retired ~2023; institutional use still possible |

**Default recommendation**: **ENA mirror** for off-cloud, **STRIDES (AWS/GCP)** for in-cloud analysis. SRA-direct only when neither is available or when SRA format itself is needed (e.g. for re-extraction of technical reads).

## Controlled-access (dbGaP) data -- check before every human accession

Everything above assumes **public** SRA data. SRA also hosts **dbGaP-controlled** human data (identifiable genomic data under a Data Use Agreement) directly alongside public runs, with the same accession format (SRR/ERR/DRR) -- nothing in the accession string itself distinguishes the two.

- **Recognize it:** the BioProject/study record on the SRA/dbGaP website shows an access column of "controlled access" (vs "public"); `prefetch`/ENA calls against a real controlled-access accession fail with an authorization error (SRA-direct) or ENA simply omits the `fastq_ftp` field for that run -- **do not treat either as a transient network failure and retry**, and do not treat a bare 403/permission error on a human accession as evidence the accession is merely mistyped.
- **This Skill is unauthenticated and out of scope for controlled-access data.** There is no public, unauthenticated path to a dbGaP-gated run. The authorized path is NCBI's `--ngc <repository-key>.ngc` flag on `prefetch` (a repository key issued per-project, per-approved-user by dbGaP after DUA approval) -- this Skill does not manage dbGaP authorization, keys, or DUAs, and an agent cannot obtain or validate one on its own.
- **The boundary:** before attempting a download of a human accession that is not already known to be from a public study (e.g. a named public cell line, a published public cohort), stop and ask the user whether the accession is dbGaP-controlled and, if so, whether they hold an approved `.ngc` repository key. Never attempt to fabricate, bypass, or work around a controlled-access restriction, and never proceed with prefetch/ENA calls on a human accession "the same way" as a public non-human one without that check.

## SRA accession hierarchy

| Prefix | Type | Granularity |
|---|---|---|
| SRR / ERR / DRR | Run | One sequencing run (file-level) |
| SRX / ERX / DRX | Experiment | Library prep + sequencing strategy |
| SRS / ERS / DRS | Sample | Biological sample |
| SRP / ERP / DRP | Study | Project (deprecated; superseded by BioProject) |
| PRJNA / PRJEB / PRJDB | BioProject | Top-level project ID |
| SAMN / SAMEA / SAMD | BioSample | Biological sample (cross-archive) |

Conversion is via SRA metadata: `pysradb metadata <ID>` or `efetch -db sra -id <UID> -format runinfo`.

The actual download unit is SRR/ERR/DRR (runs). The BioProject (PRJNA...) is the convenient top-level handle for "pull all data for paper X".

## fasterq-dump vs fastq-dump

`fasterq-dump` (sra-tools 2.10+) is the multi-threaded successor. **Always prefer it**, with two exceptions noted below.

| Aspect | fasterq-dump | fastq-dump |
|---|---|---|
| Threads | Multi (`-e N`) | Single |
| Speed | ~5-10x faster | Baseline |
| Disk overhead | Writes uncompressed FASTQ to scratch (~3x final size) | In-place; lower scratch |
| Compression | NOT built-in (post-process with pigz) | `--gzip` flag built-in |
| Single-cell technical reads | `--include-technical` works | Some 10x records need fastq-dump for full extraction |
| 10x split semantics | Sometimes incomplete | Sometimes the only way to get all reads |

The **uncompressed-scratch trap**: `fasterq-dump` writes uncompressed FASTQ first, then leaves it uncompressed. A 100 GB compressed FASTQ needs ~300 GB of scratch space + 300 GB of final output. Either compress post-hoc with `pigz` or use `--mem` to control RAM/disk tradeoff, or stick with `fastq-dump --gzip` (slower, lower scratch). Symptom of getting it wrong: "out of disk space" mid-extraction.

**No official `pigz` build exists for Windows.** On a bare Windows shell (no WSL/conda), `pigz -p N file` fails outright even though the download/extraction already succeeded. All post-hoc compression in this Skill's code patterns uses `command -v pigz` to fall back to `gzip` (single-threaded, slower, but present via Git-for-Windows/WSL/conda alike):
```bash
if command -v pigz >/dev/null 2>&1; then pigz -p "${THREADS}" file1.fastq file2.fastq
else gzip file1.fastq file2.fastq; fi
```

## prefetch and the `--max-size` trap

`prefetch` downloads `.sra` files to the configured cache before extraction. Default `--max-size 20G` silently skips runs larger than 20 GB.

```bash
# Wrong: silently skips runs >20 GB
prefetch SRR12345678

# Right: set max-size explicitly to your largest expected size
prefetch SRR12345678 --max-size 100G -p
```

For unknown-size queues, set max-size to a generous upper bound (e.g. `--max-size 200G`) or check the size first -- **but check the right size for the path you're using; `.sra` size and FASTQ size are not the same number:**

- For **`prefetch`'s `--max-size`** (governs the `.sra` download, this path): `pysradb`'s `sra_metadata(detailed=True)` **does** carry a usable field, `total_size` (checked pysradb 2.5.1) -- on real accession ERR10419835 it reported 1,744,967 bytes, matching a live `prefetch` run's actual verified `.sra` size (1,742,656 bytes) to within 0.1%. It tracks `.sra` size, not FASTQ size.
- For the **ENA mirror** (downloads FASTQ directly, a different, usually larger, file): use the ENA portal API's `fastq_bytes` field instead -- confirmed live, returns real per-mate byte counts:
```bash
curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${SRR}&result=read_run&fields=fastq_bytes&format=tsv" | tail -1
```
Sum the semicolon-separated per-file byte counts and compare against your `--max-size` bound before calling `prefetch`, or `total_size` from pysradb if going through the SRA-toolkit path.

**Before starting any download above ~10 GB** (checked via `fastq_bytes` above, or via `run_total_bases`), state the expected size and get explicit user confirmation before calling `prefetch`/`fasterq-dump`/the ENA curl loop -- this applies on top of, not instead of, the STRIDES cross-region egress-cost warning below.

## ENA mirror: direct FASTQ URLs

ENA stores FASTQ files directly (no SRA-format intermediate). Discover URLs via the ENA portal API:

```bash
curl 'https://www.ebi.ac.uk/ena/portal/api/filereport?accession=SRR12345678&result=read_run&fields=fastq_ftp,fastq_md5,read_count&format=tsv'
```

Returns TSV with semicolon-separated paired-end URLs and md5 checksums.

Direct download:
```bash
curl -O 'https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR123/078/SRR12345678/SRR12345678_1.fastq.gz'
```

ENA's mirror is typically faster than SRA's because (a) it's hosted on Aspera-aware servers, (b) the FASTQ is pre-compressed (no SRA->FASTQ conversion needed), (c) EMBL-EBI's bandwidth is generous. For most downloads in 2026, ENA is the right default.

**ENA and SRA-direct are not byte-identical, and read counts can differ.** ENA independently re-calls FASTQ from resubmitted BAM/CRAM rather than mirroring NCBI's SRA-format bytes. Checked on real accession ERR10419835: SRA-toolkit extraction (`fasterq-dump`) gives 14,052 read pairs, matching NCBI's own runinfo and `vdb-dump --info`; the ENA mirror FASTQ for the same accession has 13,911 pairs (141 fewer, ~1%), matching ENA's own `read_count` metadata exactly -- so it is not a truncated or corrupt download, the two archives just disagree. Do not validate a downloaded FASTQ's read count against a count pulled from the *other* source; validate ENA downloads against ENA's own `fastq_md5`/`read_count`, and SRA-direct downloads against `vdb-validate`/NCBI's runinfo.

## Single-cell / 10x quirks

10x Genomics records include "technical reads" (cell barcodes, UMIs) interleaved with biological reads. Default `fasterq-dump` (or `fastq-dump`) skips them, so only the cDNA file (R2) appears and CellRanger / STARsolo error out. To get all reads:

```bash
# fasterq-dump with technical reads
fasterq-dump SRR12345678 --include-technical --split-files -p -O ./fastq/

# Some 10x records require fastq-dump -- check sra-stat first
sra-stat --xml SRR12345678 | grep -E '(spotCount|baseCount|tag)'
```

For 10x v3, expect 3 files per run: R1 (barcode+UMI), R2 (cDNA), I1 (index). For 10x v2: R1 (barcode), R2 (UMI+cDNA), I1.

## MD5 / vdb-validate

Always verify downloads.

```bash
# vdb-validate for SRA-format files (toolkit path)
vdb-validate SRR12345678

# md5sum for ENA FASTQ files
md5sum -c <(echo "<expected_md5>  SRR12345678_1.fastq.gz")
```

ENA provides md5 in the portal API response. SRA-toolkit's `vdb-validate` is the equivalent for `.sra` files (different file format).

## Cloud (STRIDES) access

NCBI's STRIDES initiative mirrored SRA data to AWS Open Data (us-east-1) and GCP (us-central1). Same-region pulls have zero egress cost.

```bash
# List SRA cloud-hosted files (no NCBI auth needed)
aws s3 ls s3://sra-pub-run-odp/sra/SRR12345678/ --no-sign-request

# Direct copy to EC2 in us-east-1. The STRIDES object is named without a `.sra`
# suffix (just SRR12345678); rename on copy to keep fasterq-dump happy.
aws s3 cp s3://sra-pub-run-odp/sra/SRR12345678/SRR12345678 ./SRR12345678.sra --no-sign-request

# Then fasterq-dump locally
fasterq-dump ./SRR12345678.sra -p -e 8
```

For cloud-native analysis pipelines (Nextflow on AWS Batch, Cromwell, etc.), STRIDES is the right path.

## Code patterns

### Single SRR via ENA mirror (preferred default)

**Goal:** Download paired-end FASTQ for one SRR; verify md5; minimal dependencies.

**Approach:** Query ENA portal API for FASTQ URLs and md5; download with curl; verify with md5sum.

`echo ERR10419835 > acc.txt; bash examples/download_batch.sh acc.txt ./fastq` -- a one-line accessions file is the single-run case (ENA portal API 2.0+, curl; checked 2026-09-17). It locates the `fastq_ftp` / `fastq_md5` columns by header name (ENA always prepends `run_accession`, so a fixed `cut -f1/-f2` grabs the accession and the URL), names the missing field when one is absent (a missing `fastq_ftp` may be controlled-access, see the dbGaP section), md5-verifies every file, and lists failed accessions in `<out_dir>/failed.txt` (the script itself exits 0, so check that file).

### prefetch + fasterq-dump (SRA toolkit, classic)

`bash examples/download_single.sh <SRR> [out_dir] [threads] [max_size]` -- accepts only uppercase run accessions (`SRR`, `ERR`, or `DRR` followed by digits), then runs prefetch with explicit `--max-size`, `vdb-validate`, `fasterq-dump --split-files`, and pigz/gzip. It passes `--skip-technical`; drop that flag for 10x or other single-cell data. Toolkit cache and temporary FASTQ staging are owned directories under `out_dir` and are cleaned on every exit. If the current SRA Toolkit cannot complete a public run (including an unresolved normalized-reference dependency), the script retries the MD5-verified ENA mirror route; publication uses no-clobber hard links, so an existing file or dangling symlink is refused. The ENA helper respects the controlled-access guard: an accession without public ENA FASTQ fields produces no files and exits nonzero.

### Batch via pysradb metadata

**Goal:** Convert a list of GSE / BioProject / SRX IDs to SRR run accessions.

**Approach:** pysradb metadata returns a full hierarchy table; pull SRR column.

`python scripts/pysradb_resolve.py GSE123456 [PRJNA... SRX...]` prints the SRR runs for each ID (pysradb 2.2+, checked 2.5.1). Import `gse_to_srr` (GSE -> SRP -> SRR), `bioproject_to_runs` (full detailed metadata table) or `batch_resolve` (many IDs, failures printed and skipped) for use in a pipeline.

### Cloud (STRIDES) and 10x single-cell

`bash examples/prefetch_large.sh <SRR> [out_dir] [threads] [yes|no]` -- accepts only uppercase run accessions, checks the AWS Open Data bucket (run from EC2 in us-east-1 for zero egress), falls back to `prefetch --max-size 200G`, validates, runs `fasterq-dump --split-files`, compresses, and cleans its owned cache/staging directories. If neither Toolkit/STRIDES path completes for a public run, it safely retries the MD5-verified ENA mirror route with no-clobber publication. Pass `yes` as the 4th argument for 10x records: it swaps `--skip-technical` for `--include-technical` (10x v3 expects R1 28-bp barcode+UMI, R2 cDNA, I1 sample index).

## Common errors

| Error / symptom | Cause | Solution |
|---|---|---|
| "item not found" | Invalid accession or not in current SRA | Verify; check ENA mirror |
| Scratch disk full mid-extraction | fasterq-dump uncompressed write | Use larger scratch, `--mem`, or fastq-dump --gzip |
| Slow SRA-direct download (kbps, 9 AM-5 PM ET weekdays) | NCBI bandwidth contention; institutional users have priority | ENA or STRIDES; or run outside US business hours |
| 10x reads missing (only R2 appears) | --include-technical not set | Add the flag; verify with `sra-stat --xml` first |
| Container loses cache config; cache fills the container's small layer and the download fails | `~/.ncbi/user-settings.mkfg` is per-user, per-home and lost on rebuild | Mount a host volume at `~/.ncbi/` and persist `user-settings.mkfg`; or set `--temp` and `-O` explicitly |
| prefetch returns "success" but no file; vdb-validate / fasterq-dump then fail | Run larger than the default 20 GB `--max-size` (silent skip) | Set `--max-size` explicitly to a generous upper bound (e.g. 200G) |
| AWS bill on STRIDES | Cross-region pull (e.g. EC2 in us-west-2 against the us-east-1 bucket) | Match compute region to bucket region (us-east-1 for AWS, us-central1 for GCP) |
| Connection refused / auth fails on `ascp` against `anonftp@ftp.ncbi.nlm.nih.gov` | NCBI retired public Aspera in 2019, ENA ~2023 | Switch to HTTPS; for the fastest cloud transfer use STRIDES |
| `pigz: command not found` after a successful download | No Windows build of pigz | Fall back to `gzip` (see fasterq-dump vs fastq-dump section) |
| `curl: (6) Could not resolve host: <accession>` on the ENA path | Wrong filereport column selected (see ENA mirror code pattern) -- not a real DNS/network issue | Use header-based column lookup, not a fixed `cut -f` index |
| Authorization error / permission denied, or ENA omits `fastq_ftp`, on a human accession | Controlled-access (dbGaP) accession | Stop; see "Controlled-access (dbGaP) data" -- do not retry as a network issue |
| "Invalid run accession" | Input was not an uppercase `SRR`/`ERR`/`DRR` followed by digits | Correct the run ID; path separators, `.`, `..`, and option-like strings are rejected before download or cleanup |

## References

- NCBI. SRA Toolkit documentation. https://github.com/ncbi/sra-tools/wiki
- NCBI. STRIDES program. https://datascience.nih.gov/strides
- Leinonen R, Sugawara H, Shumway M; International Nucleotide Sequence Database Collaboration. (2011) The sequence read archive. *Nucleic Acids Res* 39:D19-D21.
- Cochrane G, Karsch-Mizrachi I, Takagi T; International Nucleotide Sequence Database Collaboration. (2016) The International Nucleotide Sequence Database Collaboration. *Nucleic Acids Res* 44:D48-D50.
- Choudhary S. (2019) pysradb: A Python package to query next-generation sequencing metadata and data from NCBI Sequence Read Archive. *F1000Research* 8:532.

## Related Skills

- entrez-search - Search the SRA db for accessions before downloading
- geo-data - GEO Series often link to SRA; gds -> sra ELink
- read-qc/quality-reports - QC the downloaded FASTQ
- read-qc/fastp-workflow - Adapter trim downloaded FASTQ
- read-alignment/bwa-alignment - Align downloaded reads
- ncbi-datasets-cli - Modern bulk path for genome data (NOT for SRA reads)
