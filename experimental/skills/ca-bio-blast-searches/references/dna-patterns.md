# DNA search patterns: blastn, megablast, dc-megablast

Code patterns for nucleotide queries against nucleotide databases. Read when writing a `blastn`/`megablast`/`dc-megablast` call. Program choice and word size: `SKILL.md`.

### Standard remote BLASTN with reproducible parameters

**Goal:** Run BLASTN with explicit, paper-quality parameters.

**Approach:** Specify program, database (refseq_select for stability), word size, expect, and a large hitlist_size to dodge the max_target_seqs trap.

**Reference (BioPython 1.83+):**
```python
from Bio.Blast import NCBIWWW, NCBIXML

handle = NCBIWWW.qblast(
    program='blastn',
    database='refseq_select_rna',
    sequence=query_seq,
    expect=1e-10,
    word_size=11,
    hitlist_size=500,  # large; filter top-N downstream
    format_type='XML',
)
record = NCBIXML.read(handle); handle.close()
top10 = sorted(record.alignments, key=lambda a: a.hsps[0].expect)[:10]
```

### Requesting megablast / dc-megablast (qblast() flag, not a program value)

**Goal:** Run high-identity DNA search (megablast) or sensitive discontiguous cross-species search (dc-megablast).

**Approach:** `NCBIWWW.qblast()` has no `program='megablast'`. Both are requested as `program='blastn'` plus a keyword flag -- confirmed live, `program='megablast'` fails immediately with `ValueError`.

**Reference (BioPython 1.83+):**
```python
# megablast (word=28, high-identity DNA, e.g. contamination screening)
handle = NCBIWWW.qblast(
    program='blastn',
    megablast=True,
    database='refseq_select_rna',
    sequence=query_seq,
    hitlist_size=500,
)

# dc-megablast (discontiguous, sensitive cross-species mRNA)
handle = NCBIWWW.qblast(
    program='blastn',
    megablast=True,
    template_type='coding',   # or 'optimal'
    template_length=18,       # 16, 18, or 21
    database='refseq_select_rna',
    sequence=query_seq,
    hitlist_size=500,
)
```
Verified live: `program='blastn', megablast=True` on a human/mouse/rat cross-species mRNA query returned 5 alignments (human + mouse only) vs. 11 from plain `blastn`/word=11 on the identical query -- reduced but non-zero cross-species sensitivity, matching the mechanism in Failure Modes below.
