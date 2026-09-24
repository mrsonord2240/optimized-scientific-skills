# scQuint: plate-based junction-cluster splicing

Use scQuint for plate-based, annotation-free junction clusters, not 10X 3'/5' data. It tests introns sharing a splice site; it does not estimate conventional cassette-exon PSI.

## Contig contract

In scQuint 0.3.3, `add_gene_annotation` prepends `chr` to GTF contigs. Therefore the supported input pair is **chr-prefixed junctions** (for example `chr1`) with an **unprefixed Ensembl-style GTF** (for example `1`). A self-consistent `1`/`1` or `chr1`/`chr1` pair can annotate zero introns. Inspect both inputs before analysis and fail fast after annotation.

```python
import numpy as np
from scquint.data import load_adata_from_starsolo, add_gene_annotation, group_introns
from scquint.differential_splicing import run_differential_splicing

adata = load_adata_from_starsolo('Solo.out/SJ/raw')
adata = add_gene_annotation(adata, 'annotation.gtf.gz')
if adata.n_vars == 0 or not adata.var['gene_id'].notna().any():
    raise ValueError('scQuint annotated zero introns: junctions must use chr-prefixed contigs and the GTF must use unprefixed contigs')
adata = group_introns(adata, by='three_prime')
adata.obs['cell_type'] = cell_types
a = np.where(adata.obs.cell_type == 'neuron')[0]
b = np.where(adata.obs.cell_type == 'glia')[0]
intron_groups, introns = run_differential_splicing(
    adata, a, b, min_cells_per_intron_group=10, min_total_cells_per_intron=10)
hits = intron_groups[intron_groups.p_value_adj < 0.05]
```

The defaults are 30 for both thresholds. With sparse data, scQuint can filter every group and raise `ValueError: not enough values to unpack` rather than return two empty tables. Treat that as insufficient support; lower both thresholds only for adequately covered plate data and report the filter change.
