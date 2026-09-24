---
name: bio-batch-downloads
description: Download large datasets from NCBI efficiently using EPost, history server, batching, rate limiting, and retry logic. Use when bulk-fetching tens of thousands of sequences, pulling all results of a large ESearch, designing reproducible pipelines, comparing E-utilities to NCBI Datasets v2 CLI, or implementing checksum-validated downloads. Encodes WebEnv TTL (~8h), EPost 200-ID limit, retmax caps, parallelization design, and integrity verification.
tool_type: python
primary_tool: Bio.Entrez
license: MIT
author: GPTomics
---

## Version Compatibility

Reference examples tested with: BioPython 1.88, NCBI Datasets CLI 18.37.0, Entrez Direct 26.0 (every runnable block verified live 2026-09-22)

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show biopython` then `help(Bio.Entrez.efetch)` to check signatures
- CLI: `datasets --version` and `efetch -version`

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Batch Downloads

**"Download N thousand records from NCBI without getting blocked"** -> The right answer is rarely "parallelize requests". For >5000 records the answer is the **history server**: search once, fetch in chunks server-side. For >100,000 records or whole genomes, the modern answer is **NCBI Datasets v2 CLI** -- the E-utilities are not optimized for bulk genome/gene data anymore.

This skill encodes (a) when to use each retrieval strategy, (b) the precise rate-limit math, (c) WebEnv lifecycle for long-running jobs, (d) how to design retry/resume, and (e) when to defect to Datasets CLI instead.

- Python: `Entrez.esearch(usehistory='y')` + chunked `Entrez.efetch()` (BioPython)
- CLI: `datasets download genome accession ...` (NCBI Datasets v2 -- preferred for genome/gene bulk)
- CLI: `epost -db <db> -input ids.txt | efetch -format fasta` (Entrez Direct)

## Required Setup

```bash
pip install biopython
# For the modern bulk path on genome/gene data:
conda install -c conda-forge ncbi-datasets-cli
```

```python
from Bio import Entrez
import time
Entrez.email = 'researcher@institution.edu'
Entrez.api_key = 'YOUR_KEY'  # 3 -> 10 req/sec; mandatory for bulk
Entrez.tool = 'project-name'
```

## Decision matrix: which retrieval strategy?

| Record count | Source | Strategy | Why |
|---|---|---|---|
| < 200 known IDs | Any db | EFetch with comma-joined `id=` | Single round-trip; trivial |
| 200-5,000 known IDs | Any db | EPost (chunked at 200) -> history -> chunked EFetch | URL length limit + chunked retrieval |
| 5,000-100,000 from a query | Any db | ESearch with `usehistory='y'` -> chunked EFetch | Push to server once; pull in batches |
| > 100,000 sequences | nucleotide/protein | Consider FTP mirror or Datasets CLI; if E-utils is still the chosen path, chunk via the history server (single stream) | NCBI throttles bulk; offline mirror is faster |
| > 1,000,000 sequences, or a literal "entire database" request | nucleotide/protein | **Refuse the EFetch-loop approach outright.** Point only to NCBI's bulk FTP/BLAST-db mirrors (`ftp.ncbi.nlm.nih.gov`) or Datasets CLI; do not attempt a chunked E-utilities loop at this scale | At this scale a single-stream E-utilities loop is a multi-week job this Skill's tools were never designed to serve -- not a "consider," a hard stop, the same way >4 concurrent workers is a hard stop |
| Whole genome assemblies | Assembly/Datasets | `datasets download genome accession ...` | Datasets v2 is the modern bulk endpoint |
| All RefSeq for a species | Datasets | `datasets download genome taxon ...` | Replaces assembly_summary.txt scraping |
| All gene records for a list | Datasets | `datasets download gene gene-id ...` | Cleaner output than EFetch gene XML |
| Raw sequencing reads | SRA | `prefetch` + `fasterq-dump` (or ENA mirror) | See `sra-data` skill |

The Datasets CLI is the right answer for any genome- or gene-centric bulk workflow as of 2023+. The E-utilities remain right for PubMed, ESummary metadata, custom queries, and anything not in the Datasets API. See `ncbi-datasets-cli` skill.

Gene-related field tags (`[GENE]`, `[Gene Name]`) require the database's official gene symbol, not a descriptive term -- e.g. `insulin[Gene Name]` returns Count=0, `INS[GENE]` does not. Verify the symbol via ESearch or the `entrez-search` skill's EInfo pattern before building a batch query around it.

## Rate-limit math (precise)

| Auth | req/sec | Sleep between calls | Bulk-friendly notes |
|---|---|---|---|
| Email only | 3 | 0.34 s | Single-threaded only; parallelism violates ToS |
| Email + API key | 10 | 0.10 s | Modest parallelism (max ~4 workers) safe |
| Institutional bulk | Negotiated | Email `eutilities@ncbi.nlm.nih.gov` | For >100K queries; courtesy expected |

NCBI's terms ask that heavy automated downloads run **outside US weekday business hours (9 AM-5 PM ET)**. Cron the job for nights/weekends; pipelines that ignore this get IP-throttled.

**Critical**: parallelizing API calls is the WRONG bulk strategy. One stream with history server + larger batches is faster AND more polite than N parallel streams. The bottleneck is rarely NCBI's throughput at small N -- it's the round-trip count.

## History server lifecycle (the long-running-job trap)

| Property | Value | Failure mode |
|---|---|---|
| TTL | 8 hours absolute (per NCBI E-utils help) | Job started Friday evening dies Saturday morning |
| Idle eviction | ~15 min empirically under load | A worker that stalls loses its WebEnv |
| Per-session isolation | One WebEnv string per session | Don't share across processes if isolation matters |
| Expired session behavior | HTTP 200 with `<ERROR>WebEnv not found</ERROR>` | Won't surface as HTTP error -- must parse body |
| Recovery | Re-run ESearch; resume at `retstart` | Need to checkpoint progress to disk |

Production pattern: checkpoint the `retstart` cursor after each successful chunk to disk; on restart, re-run ESearch (cheap), pick up `retstart` from checkpoint, continue.

## EPost specifics

EPost pushes a list of UIDs to the history server so downstream EFetch can pull by WebEnv/QueryKey instead of by ID. Two constraints:
- **200 IDs per EPost call** is the hard limit.
- **Chained posts share a WebEnv**: pass the WebEnv from the first call into subsequent calls to accumulate IDs under one session; a new QueryKey is issued per call.

To intersect: `term=#{key1} AND #{key2}` against the WebEnv produces a new key.

## Batch size guidelines per rettype

| Database | rettype | Optimal batch | Per-record payload |
|---|---|---|---|
| nucleotide | fasta | 500-1000 | ~1 KB |
| nucleotide | gb | 100-200 | ~10-50 KB |
| protein | fasta | 500-1000 | ~0.5 KB |
| protein | gp | 100-200 | ~5-30 KB |
| pubmed | medline | 1000-2000 | ~2 KB |
| pubmed | xml | 200-500 | ~10-30 KB |
| any | esummary (docsum) | 500 per call | ~1 KB |

Smaller batches for GenBank/XML because per-record payload is larger; larger batches for FASTA because the per-call HTTP overhead dominates.

## Code patterns

### Production batch fetch (history server + retry + checkpoint)

**Goal:** Download all records matching a query, robust to mid-job failures and session expiry.

**Approach:** ESearch with history; checkpoint cursor to disk; on error, retry the chunk; on session expiry (HTTP 200 with `<ERROR>` body), re-run ESearch and resume from checkpoint. The chunk is never skipped: after `max_retries` failures it raises. Bytes/str bodies are both handled.

The checkpoint stores the output file's **byte length at each committed chunk boundary**, written immediately after the chunk's bytes are flushed. Resume truncates to exactly that offset. This matters because a crash can leave *more* on disk than the checkpoint vouches for -- a chunk's records can be fully written and then the process dies before the checkpoint update. The records past the boundary are complete and well-formed, so nothing about their *content* marks them as already-counted; only the recorded offset does. Truncating to the last newline instead leaves them in place, and the resumed run appends the same records again (see the duplicate-record failure mode below).

Output is opened in binary append mode, so the fetched bytes reach disk unchanged on every platform. Opening in text mode translates NCBI's LF endings to CRLF on Windows, which silently changes a 2 MB payload's bytes and breaks any checksum comparison against the FTP manifest.

**Runnable:** `examples/robust_download.py` -> `checkpointed_download(db, term, out_path, ckpt_path, rettype='fasta', batch_size=500, max_retries=5)`. The module sets a placeholder `Entrez.email` when imported, so set yours after the import:

```python
import sys; sys.path.insert(0, 'examples')
from robust_download import checkpointed_download
from Bio import Entrez; Entrez.email = 'you@institution.edu'  # Entrez.api_key = '...' if you have one
checkpointed_download('nucleotide', 'BRCA1[GENE] AND Homo sapiens[ORGN] AND biomol_mrna[PROP] AND srcdb_refseq[PROP]',
                      'brca1_mrna.fasta', 'brca1_mrna.ckpt.json')
```

### EPost large ID list, then EFetch

**Goal:** Download by a known list of 5,000 accessions without 414 URI errors.

**Approach:** EPost in 200-ID chunks; reuse WebEnv across chunks; final fetch iterates each QueryKey by its own chunk size (retstart is relative to that key, not to the whole list).

**Runnable:** `examples/batch_by_ids.py` defines `chained_epost_fetch(db, ids, out_path, rettype='fasta', batch_size=500)` for long lists and `direct_efetch(db, ids, out_path)` for <200 IDs. It is a demo script (running or importing it fetches a 255-ID demo list), so copy the functions or edit the `small_list`/`large_list` inputs at the bottom.

### Entrez Direct from the shell

**Goal:** The same post-then-fetch flow without Python, for one-off pulls.

**Approach:** `epost` prints the WebEnv/QueryKey it created and `efetch` reads them implicitly through the pipe -- there is no `-mode webenv`. `-mode` selects the **response transport** (`text`, `xml`, `asn`, `binary`, `json`) and does not take a history-server value; the **record format** is `-format` (`fasta`, `gb`, `uilist`, ...).

Passing `-mode webenv` is not an error and does not warn: efetch falls through to a UID-list request (`rettype=uilist`), exit status stays 0, and you get an empty FASTA. Verified live against EDirect 26.0 -- `-mode webenv` gave 0 records, `-format fasta` gave 3 records / 22,122 bytes for the same input.

```bash
# Accessions or UIDs, one per line
printf 'NM_007294.4\nNM_000059.4\nNM_000546.6\n' > ids.txt

# Post once; fetch sequences through the WebEnv the pipe carries
epost -db nucleotide -input ids.txt | efetch -format fasta > out.fasta

# Verify -- 0 here means the flag spelling is wrong, not that the query was empty
grep -c '^>' out.fasta
```

For lists in the tens of thousands, stay with the Python path above: its chunking, checkpointing and expiry recovery are explicit, while EDirect's are internal.

### Integrity check after download

**Goal:** Confirm downloaded FASTA has the expected record count and no truncation.

**Approach:** Count expected (from ESearch Count) vs observed (from `SeqIO.parse`): `verify_count(path, expected)` in `examples/batch_by_ids.py` (copy it; same demo-script caveat). Skip the check when ESearch Count is 0 (no output file is written; see `examples/batch_fasta.py`).

For genome assemblies and known-checksum files, NCBI provides MD5 manifests (e.g. `md5checksums.txt` in FTP genome directories). NCBI Datasets CLI verifies checksums automatically; the FTP-direct route needs explicit `md5sum -c`.

### Compare E-utils to Datasets CLI cost

```python
def estimate_efetch_calls(total, batch_size):
    return -(-total // batch_size)  # ceiling division
```

For 100,000 nucleotide records at 500/batch with API key: 200 calls * 0.1s = 20s minimum. For the same workflow via `datasets download gene gene-id 100000`: one CLI invocation, parallel download, automatic checksum. For genome-scale bulk, Datasets wins by an order of magnitude.

### Parallelization design (modest)

**Goal:** Pull from two independent queries concurrently without violating rate limits.

**Approach:** Async with a global semaphore that enforces the API-key-permitted rate. Max 4 concurrent workers is the polite cap.

```python
import asyncio
from asyncio import Semaphore

# Pseudo-pattern; real impl needs aiohttp + Bio.Entrez async wrappers
async def fetch_with_semaphore(sem, db, id_, rettype):
    async with sem:
        # call EFetch
        await asyncio.sleep(0.1)  # rate gate
        # ... actual call

sem = Semaphore(4)
```

Never exceed 4 concurrent workers with an API key, or 1 without. Above that NCBI throttles by IP and the whole pipeline grinds.

### Regression tests

**Runnable:** `examples/test_robust_download.py` -- a network-free pytest module for the resume path and the error guards.

```bash
python -m pytest examples/test_robust_download.py -q
```

It swaps the live EFetch for a fake history server serving a locally generated record set, so it runs in ~3 s and is safe in CI. Both silent-corruption classes above are locked here: a crash leaving complete records past the committed checkpoint boundary (which produced 419 records for a 368-record query), and text-mode writes translating LF to CRLF. The fixture reproduces NCBI's blank-line-separated FASTA layout exactly, because the byte-offset arithmetic only holds if that separator is preserved.

## Failure modes

### Session expires mid-pipeline
- **Trigger:** Job runs >8h or worker idles >15 min.
- **Mechanism:** WebEnv evicted; EFetch returns HTTP 200 with `<ERROR>WebEnv not found</ERROR>` body.
- **Symptom:** Silently truncated output mid-file; downstream parsing fails on empty chunks.
- **Fix:** Parse body for `<ERROR>`; re-run ESearch and resume at checkpointed `retstart`.

### URL too long on >200 IDs
- **Trigger:** Comma-joined `id=` to EFetch with 250+ IDs.
- **Mechanism:** GET URL exceeds NCBI's ~2000 char limit.
- **Symptom:** HTTP 414 URI Too Long, or silent truncation.
- **Fix:** EPost in chunks of 200 first, then EFetch by WebEnv/QueryKey.

### Rate-limit cascade
- **Trigger:** Parallelizing without API key; or >10 req/s with key.
- **Mechanism:** NCBI returns 429; aggressive retry triggers IP-level throttle.
- **Symptom:** Pipeline gets slower and eventually stops.
- **Fix:** Add jittered exponential backoff; reduce concurrency; reach out for institutional access if bulk is the norm.

### Datasets / E-utils confusion
- **Trigger:** Building a custom assembly_summary.txt scraper instead of using Datasets.
- **Mechanism:** Datasets API is the official, supported bulk endpoint for genome/gene data; E-utils is not optimized for it.
- **Symptom:** Slow downloads, stale snapshots, missing fields.
- **Fix:** Use `datasets download genome ...` for genomes; `datasets download gene ...` for gene records. See `ncbi-datasets-cli`.

### Silent retmax cap
- **Trigger:** ESearch without `usehistory='y'`; Count > 9999.
- **Mechanism:** Legacy esearch enforces 9999 cap; the rest of the result set is silently dropped.
- **Symptom:** Batch loop terminates early; missing thousands of records.
- **Fix:** Always set `usehistory='y'` for any query expected to return >5000.

### Resume duplicates records (silent, well-formed output)
- **Trigger:** Job crashes after a chunk's bytes reach disk but before the checkpoint update lands.
- **Mechanism:** The checkpoint still names the *previous* chunk boundary, so on resume the new run re-fetches from there and appends records that are already in the file. The stale tail is complete FASTA -- last-newline truncation cannot see it, because there is no partial record to detect.
- **Symptom:** The file parses cleanly and holds **more** records than the query returned, with a duplicated run in the middle. Nothing about the file looks broken; only a count check catches it.
- **Fix:** Record the output's byte length in the checkpoint alongside the cursor, and truncate to exactly that offset on resume. Flush the chunk before writing the checkpoint, so the offset can only ever under-report -- never claim bytes that were not written. `examples/robust_download.py` does this; `examples/test_robust_download.py` locks it.

### Partial record at the end of the file
- **Trigger:** Crash mid-write, inside a chunk.
- **Mechanism:** Output file has half a record at the end.
- **Symptom:** SeqIO.parse fails on the partial record.
- **Fix:** Same byte-offset truncation as above; it discards the torn record too. Truncating to the last newline is only a fallback for checkpoints predating the byte offset, and it does not cover the duplicate case.

### Output bytes differ from the source payload
- **Trigger:** Writing the fetched body in text mode (the default `open(path, 'w')`).
- **Mechanism:** Windows translates `\n` to `\r\n` on write, so NCBI's LF-ending payload lands on disk changed. A 2.27 MB BRCA1 FASTA gained 31,994 CR bytes this way.
- **Symptom:** `md5sum` does not match the FTP manifest, or a checksum comparison fails, with no other visible difference. `SeqIO.parse` still succeeds, so the file looks fine.
- **Fix:** Open outputs with `newline=''` or in binary append mode. Verify with a byte-level check (`file`, or compare against the manifest checksum), not with `SeqIO.parse`.

## Common errors

| Error / symptom | Cause | Solution |
|---|---|---|
| HTTPError 429 | Rate limit | Sleep with backoff; get API key |
| HTTPError 414 | URL too long | EPost first |
| `<ERROR>WebEnv not found</ERROR>` (HTTP 200) | Session expired | Re-run ESearch; resume at checkpoint |
| Output file ends mid-record | Crash mid-chunk | Truncate to the checkpointed byte offset on resume |
| File parses but has extra records | Crash between chunk write and checkpoint write | Same -- the recorded offset removes the stale tail |
| Input file empty or too short | No `-id`/`-input` given | efetch reads IDs from stdin and **hangs** on a terminal; pass `-id` or pipe input, and redirect `< /dev/null` |
| Output FASTA empty, exit status 0 | `efetch -mode webenv` (invalid mode) | Use `-format fasta`; `-mode` is transport, not record format |
| Checksum mismatch, file looks fine | Text-mode write added CRLF | Open with `newline=''` or binary |
| Slow despite API key | Too few records per call | Increase batch_size to 500+ for FASTA |
| Datasets CLI faster than EFetch | Workflow is genome/gene bulk | Switch to `ncbi-datasets-cli` |

## References

- Sayers EW et al. (2024) Database resources of the National Center for Biotechnology Information in 2024. *Nucleic Acids Res* 52:D33-D43.
- Kans J. (2024) Entrez Direct: E-utilities on the Unix Command Line. NCBI Bookshelf NBK179288.
- NCBI. EPost help and Usage Guidelines. NBK25499.
- NCBI Datasets documentation: https://www.ncbi.nlm.nih.gov/datasets/docs/v2/

## Related Skills

- entrez-search - Build the query that batch-downloads will fetch
- entrez-fetch - Single-record EFetch and ESummary
- entrez-link - Chain ELink with neighbor_history for cross-db bulk
- ncbi-datasets-cli - Modern bulk endpoint for genome/gene data; preferred over E-utils for that scope
- sra-data - Raw read downloads via SRA toolkit (not via E-utilities)
- geo-data - GEO supplementary file downloads
