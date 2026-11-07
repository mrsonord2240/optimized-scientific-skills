# DepMap-style panel (Chronos)

```python
import pandas as pd, chronos
table = pd.read_csv('panel.count.txt', sep='	')       # sgRNA, Gene, then one raw-count column per sample
rep = pd.read_csv('replicatemap.txt', sep='	')        # Replicate, Sample (cell line), Control
plasmid = ['CTRL_ERS717283.plasmid']                   # plasmid-pool column(s)
days = 14                                              # days from infection to harvest
samples = list(table.columns[2:])
readcounts = table.set_index('sgRNA')[samples].T       # sequences on rows, guides on columns
readcounts.index.name = 'sequence_ID'
line = dict(zip(rep.Replicate, rep.Sample))
sequence_map = pd.DataFrame({'sequence_ID': samples})
sequence_map['cell_line_name'] = [('pDNA' if s in plasmid else line[s]) for s in samples]
sequence_map['days'] = [0 if s in plasmid else days for s in samples]
sequence_map['pDNA_batch'] = 'b1'
guide_gene_map = table[['sgRNA', 'Gene']].rename(columns={'sgRNA': 'sgrna', 'Gene': 'gene'})
neg_genes = set(pd.read_csv('NEGv1.txt', sep='	').iloc[:, 0])   # first column holds the gene symbols; or the non-targeting guides' own ids
negative_control_sgrnas = guide_gene_map[guide_gene_map.gene.isin(neg_genes)].sgrna.tolist()
model = chronos.Chronos(sequence_map={'screen': sequence_map},
                        guide_gene_map={'screen': guide_gene_map},
                        readcounts={'screen': readcounts},
                        negative_control_sgrnas={'screen': negative_control_sgrnas})
model.train(nepochs=301)
gene_effects = model.gene_effect
gene_effects.to_csv('gene_effects.csv')
```

All four inputs are dicts of DataFrames keyed by library name, not bare frames. `readcounts=`, not `reads=`; `nepochs`, not `n_steps`; `gene_effect` is an attribute. `sequence_map` needs the columns `sequence_ID`, `cell_line_name` (not `cell_line`), `days` and `pDNA_batch`; the plasmid rows are named `pDNA` with `days` 0. The constructor raises without `negative_control_sgrnas` (or `excess_variance`).

- Chronos corrects copy-number bias, screen quality and timepoints together; a single cell line does not need it.
- A separate copy-number step is `chronos.alternate_CN(gene_effect, copy_number)`, not a constructor argument.
- Install: `pip install crispr-chronos` (import name `chronos`, TensorFlow required; keep it in its own environment).

Done when `gene_effects` is written to disk.
