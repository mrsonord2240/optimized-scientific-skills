# BRIE2: Bayesian per-cell PSI

BRIE2 estimates shrunken per-cell PSI with a sequence-feature prior and tests cell covariates by likelihood-ratio test. Its PyPI sdist is broken; use the GitHub source. These install commands were not re-run for this guide.

```bash
pip install git+https://github.com/huangyh09/brie
pip install tensorflow tf_keras tensorflow-probability
# plate: sample_list.tsv is BAM path<TAB>cell ID; BAMs must be sorted and indexed
brie-count -a splicing_events.gff3 -S sample_list.tsv -o brie_counts/ -p 16
brie-quant -i brie_counts/brie_count.h5ad -c cell_metadata.tsv -o brie_quant.h5ad \
  --interceptMode gene --LRTindex All --testBase null --MCsize 3 -p 16
```

For droplet data, use `-s possorted.bam -b barcodes.tsv.gz` instead of `-S`; inspect `brie-count -h` for the installed barcode and UMI tag defaults. `briekit-event` is not a reliable event generator (`ModuleNotFoundError: parseTables`); use BRIE's precomputed exon-skipping GFF3 matching the genome build where possible.

`brie-quant` silently filters low-support events (`--minCount 50 --minUniqCount 10 --minCell 30 --minMIF 0.001` by default). Read its log before interpreting a reduced event count. With sparse but credible data, lower thresholds explicitly and report them.

```python
import scanpy as sc
import pandas as pd
a = sc.read_h5ad('brie_quant.h5ad')
lrt = pd.DataFrame({
    'ELBO_gain': a.varm['ELBO_gain'][:, 0], 'pval': a.varm['pval'][:, 0],
    'fdr': a.varm['fdr'][:, 0], 'cell_coeff': a.varm['cell_coeff'][:, 0]},
    index=a.var_names)
hits = lrt[lrt.fdr < 0.05].sort_values('fdr')
psi = a.layers['Psi']
```

Test `fdr`, not `ELBO_gain`; coefficient sign provides the covariate direction. `Psi_95CI` and `Z_std` provide uncertainty. BRIE2 does not expose a seed; repeat runs can move individual PSI values, so assess reproducibility of the significant set. For large cohorts, reduce camel-case `--batchSize` and process independent chromosome batches if memory is constrained.
