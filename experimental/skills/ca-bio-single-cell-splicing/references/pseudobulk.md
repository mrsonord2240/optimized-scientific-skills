# Replicate-aware pseudobulk junction analysis

For between-cell-type differential splicing, sum junction counts per **sample x cell type**, never into one pooled column per type. The helper below joins cells by ID, conserves reads, and rejects missing metadata or insufficient replicates.

```python
import pandas as pd

def pseudobulk_junctions(junction_counts, cell_metadata, groupby='cell_type', sample_col='sample', min_replicates=3):
    missing = junction_counts.columns.difference(cell_metadata.index)
    if len(missing):
        raise ValueError(f'{len(missing)} cells not in cell_metadata.index: {list(missing[:3])}')
    meta = cell_metadata.loc[junction_counts.columns, [groupby, sample_col]]
    if meta.isna().any().any():
        raise ValueError(f'NaN in {groupby!r} or {sample_col!r}: those cells would be dropped')
    pb = junction_counts.T.groupby([meta[groupby].to_numpy(), meta[sample_col].to_numpy()]).sum().T
    assert pb.to_numpy().sum() == junction_counts.to_numpy().sum()
    groups = pd.DataFrame({'sample': [f'{g}__{s}' for g, s in pb.columns],
                           'group': [g for g, _ in pb.columns]})
    n_rep = groups.group.value_counts()
    if (n_rep < min_replicates).any():
        raise ValueError(f'need >= {min_replicates} samples per {groupby}; got {n_rep.to_dict()}')
    pb.columns = groups['sample']
    return pb, groups
```

Leafcutter expects clustered intron counts (`*_perind_numers.counts.gz` layout), not arbitrary junction rows. For a two-group pair:

```bash
leafcutter_ds.R -i 5 -g 3 -c 20 -o pb_ds pb_counts.txt.gz pb_groups.txt
```

Defaults require support across replicates. Lowering `-g` and `-i` to one yields p-values without biological replication; do not present those as differential-splicing evidence. rMATS consumes BAMs, so merge cell BAMs per sample x cell type first.
