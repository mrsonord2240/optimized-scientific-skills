# Protein search patterns: organism restriction, short peptides

Code patterns for protein queries. Read when writing a `blastp` call that restricts by organism or searches a peptide under 30 aa. Matrix / gap-cost / CBS choices: `SKILL.md`.

### Protein search with organism restriction

**Goal:** Find mammalian homologs of a query protein in Swiss-Prot.

**Approach:** `entrez_query` filters the BLAST search space pre-execution; faster and more meaningful E-values than post-filtering.

**Reference (BioPython 1.83+):**
```python
handle = NCBIWWW.qblast(
    program='blastp',
    database='swissprot',
    sequence=protein_seq,
    entrez_query='Mammalia[Organism]',
    expect=1e-5,
    composition_based_statistics=2,
    hitlist_size=200,
)
record = NCBIXML.read(handle); handle.close()
```

### Short peptide search

```python
handle = NCBIWWW.qblast(
    program='blastp',
    database='swissprot',
    sequence=peptide_seq,  # < 30 aa
    matrix_name='PAM30',
    word_size=2,
    gapcosts='9 1',  # required: PAM30 rejects BLOSUM62's default gap costs (11,1) -- see word-size/gap-cost table above
    expect=1000,  # short queries need permissive cutoff
    composition_based_statistics=3,
    hitlist_size=100,
)
```
Verified live: omitting `gapcosts` raises `ValueError: Error message from NCBI: ... Gap existence and extension values of 11 and 1 not supported for PAM30`; with `gapcosts='9 1'` the call succeeds.
