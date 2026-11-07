# Ids



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

**Approach:** Count expected (from ESearch Count) vs observed (from `SeqIO.parse`): `verify_count(path, expected)` in `examples/batch_by_ids.py` (copy it; same demo-script caveat). Skip the check when ESearch Count is 0 and no output file is written.

For genome assemblies and known-checksum files, NCBI provides MD5 manifests (e.g. `md5checksums.txt` in FTP genome directories). NCBI Datasets CLI verifies checksums automatically; the FTP-direct route needs explicit `md5sum -c`.

### Compare E-utils to Datasets CLI cost

```python
def estimate_efetch_calls(total, batch_size):
    return -(-total // batch_size)  # ceiling division
```

For 100,000 nucleotide records at 500/batch with API key: 200 calls * 0.1s = 20s minimum. For the same workflow via `datasets download gene gene-id 100000`: one CLI invocation, parallel download, automatic checksum. For genome-scale bulk, Datasets wins by an order of magnitude.
