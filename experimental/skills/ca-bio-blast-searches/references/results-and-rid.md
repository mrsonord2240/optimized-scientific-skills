# Results handling and RID polling

Saving XML, extracting hits, and polling a long job by RID. Read when parsing results or when a search will outlive one blocking `qblast()` call. RID lifecycle table: `SKILL.md`.

### Save XML for re-parsing

`scripts/blast_rid.py` saves XML; `scripts/parse_results.py` reads one result with `NCBIXML.read()` and writes the first HSP of each alignment as TSV.

### Hit extraction with identity + coverage filtering

**Goal:** Return structured top hits with biological metrics, not just E-values.

**Approach:** Walk alignments + first HSP; compute identity and query coverage as fractions; sort by bit-score (database-size invariant) not E-value.

Filter the TSV by identity and coverage, and sort by bit score. The source's coverage calculation uses alignment length divided by query length; gap and translated-search units require special care.

### Programmatic RID polling for long jobs

**Goal:** Submit a long job, keep the RID, poll without blocking, resume later.

**Approach:** `qblast()` never exposes the RID, so use the NCBI BLAST URL API directly (`scripts/blast_rid.py`, stdlib only): `Put` returns RID + RTOE, `SearchInfo` returns `Status=WAITING|READY|FAILED|UNKNOWN`, `Get` returns the XML. It waits RTOE, polls at most once per 60 s, refuses queries without a defline, and exits with the NCBI error text on a rejected submit.
```bash
# from the Skill directory
python scripts/blast_rid.py run --query q.fa --program blastn --db refseq_select_rna --hitlist 500 --expect 1e-10 --out hits.xml
python scripts/blast_rid.py submit --query q.fa --program tblastn --db nr      # prints RID; later: status RID / fetch RID --out hits.xml
```
The saved XML parses with `NCBIXML.read()` (see `scripts/parse_results.py`).
