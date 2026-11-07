# Python only

```python
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

dds = DeseqDataSet(counts=count_df, metadata=metadata, design='~condition')   # count_df: samples x genes
dds.deseq2()
stats = DeseqStats(dds, contrast=('condition', 'treated', 'control'))         # (column, numerator, reference)
stats.summary()
stats.results_df.to_csv('results.csv')
```

- `count_df` has samples as rows and genes as columns, the transpose of the R layout.
- Batch or pairing: `design='~batch + condition'`.
- PyDESeq2 has no likelihood-ratio test. For "any difference across three or more levels", use R (`routes/multi-level.md`).
- If R with DESeq2 is installed, prefer the R routes: they are one command.

Done when `results.csv` exists and you have answered from it.

More detail: `references/pydeseq2.md`.
