# Search

taxonomy (the most-misused BLAST setting)

For `blastn`, the `-task` flag picks among heuristics with different word sizes and gap parameters.

| `-task` | Word | Gapped | Use case | Mistake to avoid |
|---|---|---|---|---|
| `megablast` (default) | 28 | linear | >=95% identity, intra-species, primer hits, contamination check | Used for cross-species and misses everything |
| `dc-megablast` | 11 (discontiguous) | yes | Cross-species mRNA homology | Underused -- this is what `blastn` "should" be for cross-species |
| `blastn` | 11 | yes | General sensitive DNA | Slower than dc-megablast for same job |
| `blastn-short` | 7 | yes | Queries <50 nt (primers, small RNAs) | Default megablast can't seed at length 7 |
| `rmblastn` | 11 | yes | Repeat masking; bundled with RepeatModeler | Specialized |

For `blastp`:

| `-task` | Word | Use case |
|---|---|---|
| `blastp` (default) | 3 | General protein similarity |
| `blastp-fast` | 6 | Faster, less sensitive |
| `blastp-short` | 2 | Peptides <30 aa, with PAM30 + word_size=2 typical |

## Soft vs hard masking

| Setting | Effect on seed | Effect on extension | Effect on score |
|---|---|---|---|
| `-soft_masking true` (default for several tasks) | Skip masked positions when seeding | Allow extension through masked | Score includes masked positions |
| `-soft_masking false` + `-dust yes` / `-seg yes` | Skip masked positions when seeding | Skip masked positions in extension | Score excludes masked positions |
| Hard-mask in input FASTA (N or X) | Hard exclusion everywhere | Hard exclusion | Treated as mismatches |

Soft masking is correct for almost all cases. Hard masking creates artificial mismatches at masked boundaries and can split true alignments. The exception: searching against a database of repeats explicitly, where hard masking on the query is the right choice.

## Thread scaling

BLAST+ parallelizes per-query (with `-num_threads`) but is I/O bound past ~16 threads on most hardware. For >100,000 query batches the better answer is splitting the input FASTA into N chunks and running N parallel `blastn` invocations -- this saturates CPUs better than `-num_threads 64`.

| Threads | Typical speedup vs single | Notes |
|---|---|---|
| 1-8 | Near-linear | Default sweet spot |
| 8-16 | Sub-linear (1.5-2x over 8) | Useful on big SMP boxes |
| 16-32 | Diminishing returns | I/O bound for most DBs |
| 32+ | Often slower | Cache thrash + I/O contention |

For massive workflows, prefer **DIAMOND** (Buchfink et al. 2021 *Nat Methods* 18:366) or **MMseqs2** (Steinegger & Soding 2017 *Nat Biotechnol* 35:1026) -- 100-10,000x faster than BLASTP at comparable sensitivity. See `bio-remote-homology` skill.

## Output format reference (`-outfmt`)

| `-outfmt` | Description | Use |
|---|---|---|
| 0 | Pairwise (default; human-readable) | Debugging, inspection |
| 5 | XML | Programmatic parsing (Bio.SearchIO) |
| 6 | Tabular (no header) | Most pipelines |
| 7 | Tabular with comment headers | Self-documenting |
| 11 | ASN.1 binary | Re-parse with later versions |

Custom tabular fields:
```bash
blastn -query q.fa -db db -outfmt "6 qseqid sseqid pident length qcovs qcovhsp evalue bitscore staxids sscinames stitle"
```

Field key fields for analysis:
- `pident` = percent identity over the HSP (NOT the query); for query-level, use `qcovhsp`
- `qcovs` = total query coverage by all HSPs of this subject (the "coverage" most users want)
- `qcovhsp` = query coverage by best HSP alone (use when there's only one HSP per hit)
- `staxids` = taxonomy IDs assigned to subject sequences; taxonomy filtering still requires v5 plus its companion taxonomy data

## Prebuilt NCBI databases via `update_blastdb.pl`

```bash
# List available
update_blastdb.pl --showall pretty | grep -E 'refseq|swissprot|nt|nr'

# Download (with decompress)
update_blastdb.pl --decompress refseq_select_rna

# Download specific volume of split database
update_blastdb.pl --decompress refseq_protein

# Download with parallelism
update_blastdb.pl --decompress --num_threads 4 refseq_select_rna
```

Sizes (approximate, 2026):
- `refseq_select_rna`: ~5 GB
- `refseq_protein`: ~30 GB
- `swissprot`: <1 GB
- `nt`: ~250 GB
- `nr`: ~300 GB

For most use cases, `refseq_select_*` is the right starting point. `nt`/`nr` are storage-heavy and reproducibility-hostile. Downloads past a few GB can take hours -- run large `update_blastdb.pl` pulls (`refseq_protein` and up) overnight rather than blocking on them interactively.

## Code patterns

### Build and search a custom protein database

**Goal:** Build a BLAST+ protein database from a custom FASTA and search against it.

**Approach:** `makeblastdb` with v5 + parse_seqids + hash_index; `blastp` with explicit outfmt.

**Reference (NCBI BLAST+ 2.15+):**
```bash
#!/bin/bash
# Reference: NCBI BLAST+ 2.15+ | Verify API if version differs

REF=reference_proteins.fasta
DB=ref_prot_db
QUERY=query.fasta
OUT=hits.tsv

makeblastdb -in "$REF" -dbtype prot \
            -blastdb_version 5 -parse_seqids -hash_index \
            -title "$REF $(date +%Y-%m-%d)" \
            -out "$DB"

blastp -query "$QUERY" -db "$DB" \
       -evalue 1e-10 \
       -num_threads 8 \
       -max_target_seqs 500 \
       -outfmt "6 qseqid sseqid pident length qcovs evalue bitscore stitle" \
       -out "$OUT"

# Top hit per query by bit-score (column 7)
sort -k1,1 -k7,7gr "$OUT" | awk '!seen[$1]++' > top_hit_per_query.tsv
```

### Cross-species DNA with dc-megablast

```bash
blastn -query mouse_cdna.fa -db human_refseq_rna \
       -task dc-megablast \
       -word_size 11 \
       -evalue 1e-10 \
       -outfmt "6 qseqid sseqid pident length qcovs evalue bitscore" \
       -num_threads 8 \
       -out cross_species.tsv
```

### Short primer search

```bash
blastn -query primers.fa -db genome_db \
       -task blastn-short \
       -word_size 7 \
       -evalue 1000 \
       -outfmt 6 \
       -out primer_hits.tsv
```

### Taxonomy-filtered search (BLAST v5 only)

```bash
# Restrict to specific taxids
blastp -query query.fa -db nr \
       -taxids 9606,10090,10116 \
       -outfmt "6 qseqid sseqid staxids sscinames evalue bitscore" \
       -out mammalian_hits.tsv

# Or to a taxid subtree (NCBI BLAST+ 2.13+)
echo 9606 > human_only.txt
blastp -query query.fa -db nr -taxidlist human_only.txt -outfmt 6 -out human_hits.tsv
```

### Extract subject sequences for top hits

```bash
# Requires database built with -parse_seqids
cut -f2 top_hit_per_query.tsv | sort -u > hit_accessions.txt
blastdbcmd -db ref_prot_db -entry_batch hit_accessions.txt -out hits.fasta

# Pull a range of a sequence
blastdbcmd -db genome_db -entry NC_000001.11 -range 1000000-1001000 -out region.fa
```

### Reciprocal best hit (RBH) for ortholog candidates

See `bio-ortholog-inference` skill for the principled treatment. Quick version:

```bash
blastp -query A.fa -db B_db -outfmt 6 -evalue 1e-5 -num_threads 8 \
       -max_target_seqs 5 -out A_vs_B.tsv
blastp -query B.fa -db A_db -outfmt 6 -evalue 1e-5 -num_threads 8 \
       -max_target_seqs 5 -out B_vs_A.tsv

# Best forward + reverse, intersect
awk '!seen[$1]++ {print $1"\t"$2}' A_vs_B.tsv | sort > A_best
awk '!seen[$1]++ {print $1"\t"$2}' B_vs_A.tsv | sort > B_best
awk 'NR==FNR{a[$1]=$2; next} a[$2]==$1' A_best B_best > rbh.tsv

# rbh.tsv columns are (B_accession, A_accession) -- column 1 came from B_best's query
# column (the B-vs-A search), column 2 from its best-hit column (an A accession).
# Extract both sides from the matching DB (requires -parse_seqids at build time):
cut -f1 rbh.tsv | sort -u > rbh_B_accessions.txt
cut -f2 rbh.tsv | sort -u > rbh_A_accessions.txt
blastdbcmd -db B_db -entry_batch rbh_B_accessions.txt -out rbh_B_hits.fasta
blastdbcmd -db A_db -entry_batch rbh_A_accessions.txt -out rbh_A_hits.fasta
```

This works but does NOT handle paralog mis-pairs from gene duplication; for that use OrthoFinder or OMA (in `bio-ortholog-inference`).

### Python wrapper with version pinning

```python
import subprocess
import shutil


def require_tool(name, min_version=None):
    if not shutil.which(name):
        raise RuntimeError(f'{name} not on PATH')
    out = subprocess.run([name, '-version'], capture_output=True, text=True)
    print(f'  {out.stdout.strip().splitlines()[0]}')


def run_blast(query, db, out, program='blastp', evalue=1e-10, threads=8, hitlist=500):
    require_tool(program)
    cmd = [program, '-query', query, '-db', db, '-out', out,
           '-evalue', str(evalue),
           '-num_threads', str(threads),
           '-max_target_seqs', str(hitlist),
           '-outfmt', '6 qseqid sseqid pident length qcovs qcovhsp evalue bitscore stitle']
    subprocess.run(cmd, check=True)


def parse_tabular(path):
    cols = ['qseqid', 'sseqid', 'pident', 'length', 'qcovs', 'qcovhsp', 'evalue', 'bitscore', 'stitle']
    rows = []
    with open(path) as f:
        for line in f:
            vals = line.rstrip('\n').split('\t')
            d = dict(zip(cols, vals))
            for k in ('pident', 'qcovs', 'qcovhsp', 'evalue', 'bitscore'):
                d[k] = float(d[k])
            d['length'] = int(d['length'])
            rows.append(d)
    return rows
```

## Failure modes

### `nt`/`nr` size shock
- **Trigger:** `update_blastdb.pl --decompress nt` without realizing the size.
- **Mechanism:** `nt` is ~250 GB compressed, ~1 TB indexed.
- **Symptom:** Disk fills mid-download; partial DB unusable.
- **Fix:** Use `refseq_select` for most workflows; only pull `nt`/`nr` with intent and >1 TB free.

### Missing `-parse_seqids`
- **Trigger:** Built DB without `-parse_seqids`; later try `blastdbcmd -entry`.
- **Mechanism:** Without the parsed index, `blastdbcmd` can't look up by accession.
- **Symptom:** `Error: ... not found in database`.
- **Fix:** Rebuild with `-parse_seqids` (cheap if FASTA still on disk).

### Wrong `-task` for the question
- **Trigger:** Default `blastn` for cross-species mRNA (word=11 but ungapped seeding).
- **Mechanism:** Discontiguous seed (`dc-megablast`) is much more sensitive across species.
- **Symptom:** Far fewer hits than the question warrants.
- **Fix:** Use `-task dc-megablast` for cross-species; `-task megablast` only for >=95% identity.

### Thread saturation
- **Trigger:** `-num_threads 64` on a 32-core box.
- **Mechanism:** I/O bound past ~16 threads; cache thrash hurts past CPU count.
- **Symptom:** No speedup or slowdown.
- **Fix:** Cap at 8-16; split FASTA and run parallel processes instead for very large batches.

### Taxonomy filtering: v5 no-op versus v4 hard error
- **Trigger:** `-taxids`/`-taxidlist` on a v5 DB before `taxdb.tar.gz` is in `$BLASTDB`/CWD, or on a v4 DB.
- **Mechanism:** v5 stores the taxid assignment and can filter, but needs `taxdb.tar.gz` for taxonomy lookup plus an auto-fetched `taxonomy4blast.sqlite3` on first use. v4 can retain assigned taxids but does not implement taxonomy filtering.
- **Symptom:** v5 without the data prints `requires additional data files` yet exits 0 and returns the full, unfiltered hit set. A v4 DB with taxonomy data present fails loudly with a nonzero exit and `Taxonomy filtering is not supported in v4 BLAST dbs`.
- **Fix:** For v5, fetch `taxdb.tar.gz` into `$BLASTDB` or the CWD and confirm filtering by row count, not exit code. For v4, rebuild as v5; fetching taxonomy data cannot make a v4 DB filterable.

### Soft-masking confusion
- **Trigger:** Hard-masking input (replacing repeats with N or X) instead of using `-dust`/`-seg`.
- **Mechanism:** Hard-mask creates artificial mismatches at boundaries.
- **Symptom:** True alignments split into multiple short HSPs.
- **Fix:** Pass unmasked FASTA + soft-mask via `-soft_masking true` + `-dust yes`/`-seg yes`.

### `max_target_seqs` truncation
- **Trigger:** `-max_target_seqs 10` (Shah et al. 2019 *Bioinformatics* 35:1613).
- **Mechanism:** Early termination, not top-N filter.
- **Symptom:** Different top-10 than `-max_target_seqs 500` + post-filter.
- **Fix:** Set `-max_target_seqs` large (500+); filter top N in awk/Python.

## Common errors

| Error / symptom | Cause | Solution |
|---|---|---|
| `BLAST Database error` | DB path wrong, or alias missing | `blastdbcmd -db <db> -info` to confirm |
| `Error: entry not found` | Built without `-parse_seqids` | Rebuild |
| Taxonomy filter no-op (exit 0, unfiltered rows) | v5 DB missing `taxdb.tar.gz` | Fetch `taxdb.tar.gz` into `$BLASTDB` or the CWD; exit 0 does not mean it filtered, so check row counts |
| `Taxonomy filtering is not supported in v4 BLAST dbs` (nonzero exit) | v4 DB cannot apply `-taxids`/`-taxidlist` | Rebuild as v5; v4 may report assigned taxids but cannot filter on them |
| Threads >16 not faster | I/O bound | Split input + parallel invocations |
| `nt` download fills disk | Database is huge | Use refseq_select |
| `Sequence too short` | Query < word_size | Use `-task blastn-short` (word=7) |
| Out of memory | Single large query | Reduce `-num_threads`, split query |
| Diffing raw `-out` file against an awk/sort-derived file shows spurious differences | Native BLAST `-out` files use CRLF on Windows; shell pipelines emit LF | Normalize with `tr -d '\r'` before diffing or scripting against a raw `-out` file |
