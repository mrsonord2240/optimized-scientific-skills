# SRA retrieval

Read when converting an SRA UID to SRR accessions and run metrics, or fetching SRA XML.

Assumes the Required Setup block in `SKILL.md` (`Entrez.email`, imports, `expect_start`).

### sra

| rettype | retmode | Returns | Use when |
|---|---|---|---|
| `runinfo` | `text` | CSV of run metadata | Convert SRA UID -> SRR accession + Run metrics |
| `xml` | `xml` | Full SRA XML hierarchy | Need BioSample/BioProject linkage in one call |

### SRA UID -> SRR accession + run metrics

**Goal:** Convert an opaque SRA UID into the SRR run accession plus Bases/Spots metrics, in one EFetch.

**Approach:** `rettype='runinfo'` returns a CSV row per run.

```python
def sra_runinfo(uids):
    h = Entrez.efetch(db='sra', id=','.join(uids), rettype='runinfo', retmode='text')
    raw = h.read(); h.close()
    # db='sra' returns bytes here despite retmode='text' (Biopython 1.88) -- decode first.
    text = raw.decode() if isinstance(raw, bytes) else raw
    lines = expect_start(text, 'Run,').strip().split('\n')
    header = lines[0].split(',')
    return [dict(zip(header, row.split(','))) for row in lines[1:]]
```
