# Database

## Database format: v5 vs v4

NCBI introduced BLAST database v5 in BLAST+ 2.10 (2020). Both v4 and v5 can store per-sequence taxid assignments supplied through `-taxid`/`-taxid_map` at build time. **Only v5 supports `-taxids`/`-taxidlist` filtering.** A v5 filter also needs NCBI's `taxdb.tar.gz` (scientific-name lookup: `taxdb.bti`/`taxdb.btd`) present in `$BLASTDB` or the working directory, and `blastp`/`blastn` auto-fetch a further ~98 MB `taxonomy4blast.sqlite3` over the network on the first `-taxids`/`-taxidlist` call. Checked on NCBI BLAST+ 2.17.0+.

| Feature | v4 | v5 |
|---|---|---|
| Default for prebuilt NCBI dbs | No (legacy) | Yes (since 2020) |
| Per-sequence taxid at build time (`-taxid_map`) | Yes (stored, but not filterable) | Yes |
| `-taxids`, `-taxidlist` support (with `taxdb.tar.gz` present) | No | Yes |
| `blastdbcmd -taxids` | No | Yes |
| New `-info` output fields | No | Yes |

`update_blastdb.pl` downloads v5 by default. When building a database manually with `makeblastdb`, v5 format requires `-blastdb_version 5`. **Always pass `-blastdb_version 5` and `-parse_seqids` when building from scratch.**

**A v5 taxonomy filter silently no-ops without `taxdb.tar.gz`.** A v5 DB built with `-taxid_map` alone is NOT enough for `-taxids`/`-taxidlist` to actually filter. Without `taxdb.tar.gz` in `$BLASTDB` or the CWD, `blastp`/`blastn` print `The -taxids command line option requires additional data files ...` to stderr but still **exit 0 and return the unfiltered hit set** -- verified on 2.17.0+: a `-taxidlist` restricted to one taxid returned both the matching and the non-matching hit, byte-identical to the unfiltered run, at exit code 0.

```bash
# Fetch once, into $BLASTDB or the working directory
curl -O https://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz
tar -xzf taxdb.tar.gz   # taxdb.bti, taxdb.btd

# The first -taxids/-taxidlist call also auto-fetches ~98 MB taxonomy4blast.sqlite3
# over the network -- expect a one-time delay, not a hang.
```

**Detect the v5 no-op -- exit 0 does not mean it worked.** Compare row counts before and after filtering, or grep stderr for `requires additional data files`. If a `-taxidlist` search returns the same row count as the unfiltered search, the filter did not apply: fetch `taxdb.tar.gz` and rerun before trusting the result as filtered. (Verified: 2 unfiltered rows -> still 2 rows without `taxdb.tar.gz` -> 1 row, correctly excluding the non-matching taxid, once `taxdb.tar.gz` was present.)

**Do not treat this as a v4 fallback.** Although a v4 database can retain and report assigned taxids (for example through `blastdbcmd -outfmt %T` or `staxids`), it cannot apply `-taxids`/`-taxidlist`. With the taxonomy data installed, BLAST+ fails loudly (nonzero exit; `Taxonomy filtering is not supported in v4 BLAST dbs`) rather than returning a silently unfiltered result. Rebuild as v5 for taxonomy filtering.

## `makeblastdb` flag taxonomy

| Flag | Effect | When |
|---|---|---|
| `-dbtype nucl` or `-dbtype prot` | Required | Always |
| `-parse_seqids` | Indexes accessions so `blastdbcmd -entry <acc>` works | Almost always (downstream extraction) |
| `-hash_index` | Speeds up extraction by accession | Large dbs |
| `-blastdb_version 5` | Use v5 format | Always |
| `-taxid 9606` | Single taxid for all seqs | Single-species DB |
| `-taxid_map file.tsv` | Per-sequence taxid mapping (seqid<TAB>taxid) | Multi-species DB |
| `-mask_data masking.asnb` | Apply precomputed soft-masking | Production pipelines |
| `-title "..."` | Free-text label | Cosmetic |
| `-out path/prefix` | DB file path prefix | Always |

```bash
makeblastdb -in reference.fasta -dbtype nucl \
            -blastdb_version 5 \
            -parse_seqids \
            -hash_index \
            -title "Custom reference 2026-05" \
            -out custom_db
```

## `-task`
