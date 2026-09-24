# Streaming Large Stockholm Databases

Read this reference for Pfam-A.full, BFD, or another alignment collection too large to load into memory.

`Bio.AlignIO.read()` is in-memory; for Pfam-A.full (multi-gigabyte; ~22,000 family alignments in Pfam 37) or BFD (>2 TB), use `pyhmmer.easel.MSAFile` for streaming Stockholm or A2M (checked on pyhmmer 0.12.3).

```python
import pyhmmer

with pyhmmer.easel.MSAFile('Pfam-A.full', digital=True) as msa_file:
    for msa in msa_file:
        nseq, alen = len(msa.sequences), len(msa.alignment[0])   # DigitalMSA has no .nseq / .alen
        if nseq < 50:
            continue
        weights = msa.compute_weights(method='pb')
        print(msa.name, nseq, alen, f'sum_w={sum(weights):.1f}')   # msa.name is str, not bytes
```

`msa.compute_weights(method='pb')` computes Henikoff PB weights via the same Easel routine HMMER uses; the weights sum to the number of sequences (not Neff). For an Henikoff-style Neff estimate, see `msa-parsing/examples/neff.py`.
