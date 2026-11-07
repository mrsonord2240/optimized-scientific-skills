# History



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
- **Fix:** Use `datasets download genome ...` for genomes; `datasets download gene ...` for gene records. See `bio-ncbi-datasets-cli`.

### Silent retmax cap
- **Trigger:** ESearch without `usehistory='y'`; Count > 9999.
- **Mechanism:** Legacy esearch enforces 9999 cap; the rest of the result set is silently dropped.
- **Symptom:** Batch loop terminates early; missing thousands of records.
- **Fix:** Always set `usehistory='y'` for any query expected to return >5000.

### Resume duplicates records (silent, well-formed output)
- **Trigger:** Job crashes after a chunk's bytes reach disk but before the checkpoint update lands.
- **Mechanism:** The checkpoint still names the *previous* chunk boundary, so on resume the new run re-fetches from there and appends records that are already in the file. The stale tail is complete FASTA -- last-newline truncation cannot see it, because there is no partial record to detect.
- **Symptom:** The file parses cleanly and holds **more** records than the query returned, with a duplicated run in the middle. Nothing about the file looks broken; only a count check catches it.
- **Fix:** Record the output's byte length in the checkpoint alongside the cursor, and truncate to exactly that offset on resume. Flush the chunk before writing the checkpoint, so the offset can only ever under-report -- never claim bytes that were not written. `scripts/robust_download.py` implements this.

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
| Datasets CLI faster than EFetch | Workflow is genome/gene bulk | Switch to `bio-ncbi-datasets-cli` |
