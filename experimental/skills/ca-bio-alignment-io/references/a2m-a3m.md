# A2M / A3M Conventions

Read this reference when working with HMMER, HH-suite, or ColabFold alignment output.

A2M (HMMER) and A3M (HHsuite, ColabFold) encode match vs insert columns by case (uppercase / `-` = match column, lowercase = insert; `.` pads insert columns in padded A2M). HMMER `hmmalign --outformat A2M` writes **ragged** rows (lowercase inserts, no `.` padding), so `AlignIO.read(..., 'fasta')` raises `ValueError: Sequences must all be the same length`. A2M padded by HH-suite loads as a rectangular MSA (`AlignIO.read(..., 'fasta')` or `Align.read(..., 'a2m')`). A3M is not padded and has no reader in AlignIO or pyhmmer (`format='a3m'` raises `InvalidParameter`; pyhmmer reads `a2m`), so convert it first with HHsuite `reformat.pl a3m a2m in.a3m out.a2m`, or take the match columns straight from the ragged rows.

**reformat.pl pitfall:** HHsuite's `reformat.pl a3m a2m` uses the FIRST sequence in the A3M as the match-state reference (verified: the same A3M with a different first record gives different padding/case for the other rows). ColabFold MSAs typically place the query first, which is the desired reference; merged or sorted A3Ms can have a non-query first sequence, producing match-state assignments that mis-align the query. Move the query record to the first position before reformatting (`hhfilter` does not reorder records, so it is not a fix). A3M files emitted by `hhblits` always have the query first; A3M files concatenated from MSA databases do not.

```python
from Bio import SeqIO

# SeqIO.parse does not require equal row lengths, so this works for ragged hmmalign A2M and padded A2M alike
match_only_seqs = {
    r.id: ''.join(c for c in str(r.seq) if c.isupper() or c == '-')
    for r in SeqIO.parse('hits.a2m', 'fasta')
}   # every row has one character per model match state (117 for the Pfam PF00042 HMM)
```
