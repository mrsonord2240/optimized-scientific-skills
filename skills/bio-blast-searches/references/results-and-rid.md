# Results handling and RID polling

Saving XML, extracting hits, and polling a long job by RID. Read when parsing results or when a search will outlive one blocking `qblast()` call. RID lifecycle table: `SKILL.md`.

### Save XML for re-parsing

`examples/save_and_parse.py` (`run_and_save()` writes `qblast()` output to disk, `parse_hits()` re-reads it with `NCBIXML.read()`).

### Hit extraction with identity + coverage filtering

**Goal:** Return structured top hits with biological metrics, not just E-values.

**Approach:** Walk alignments + first HSP; compute identity and query coverage as fractions; sort by bit-score (database-size invariant) not E-value.

`examples/basic_blast.py` `top_n_by_bitscore(record, n, min_identity, min_coverage)` (protein version with coverage-only filter: `examples/blastp_filtered.py` `filter_top`).

### Programmatic RID polling for long jobs

**Goal:** Submit a long job, keep the RID, poll without blocking, resume later.

**Approach:** `qblast()` never exposes the RID, so use the NCBI BLAST URL API directly (`scripts/blast_rid.py`, stdlib only): `Put` returns RID + RTOE, `SearchInfo` returns `Status=WAITING|READY|FAILED|UNKNOWN`, `Get` returns the XML. It waits RTOE, polls at most once per 60 s, refuses queries without a defline, and exits with the NCBI error text on a rejected submit.
```bash
# from the Skill directory
python scripts/blast_rid.py run --query q.fa --program blastn --db refseq_select_rna --hitlist 500 --expect 1e-10 --out hits.xml
python scripts/blast_rid.py submit --query q.fa --program tblastn --db nr      # prints RID; later: status RID / fetch RID --out hits.xml
```
The saved XML parses with `NCBIXML.read()` (see `examples/save_and_parse.py`).
