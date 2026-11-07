# PyDESeq2 (Python alternative)

## PyDESeq2 (Python alternative)

```python
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

dds = DeseqDataSet(counts=count_df, metadata=metadata, design='~condition')
dds.deseq2()
stat_res = DeseqStats(dds, contrast=('condition', 'treated', 'control'))
stat_res.summary()
results_df = stat_res.results_df
```

PyDESeq2 0.5+ supports Wald, multi-factor designs, and apeglm shrinkage. No LRT yet. Results are numerically close to R DESeq2 but with small differences from MLE solver choice.
