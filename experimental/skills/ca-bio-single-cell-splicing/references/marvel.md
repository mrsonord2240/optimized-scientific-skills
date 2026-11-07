# MARVEL: plate-based PSI and cell-type comparisons

Use MARVEL for full-length plate data. Build a **wide** STAR splice-junction matrix: first column `coord.intron` formatted `chr:start:end`; remaining columns are cell IDs and integer unique-junction counts. Long-format junction tables are not accepted.

```r
library(MARVEL); library(data.table); library(Seurat)

# `gtf` must be global: Preprocess_rMATS ignores its GTF= argument in MARVEL 2.0.5.
gtf <- fread('annotation.gtf', header=FALSE, sep='\t', quote='', data.table=FALSE)
for (ev in c('SE', 'MXE', 'RI', 'A5SS', 'A3SS')) {
  tab <- Preprocess_rMATS(read.table(sprintf('rmats_out/fromGTF.%s.txt', ev), header=TRUE, sep='\t'),
                           GTF=gtf, EventType=ev)
  write.table(tab, sprintf('events_%s.txt', ev), sep='\t', quote=FALSE, row.names=FALSE)
}

sj_files <- list.files('star_pass2/', pattern='SJ.out.tab$', full.names=TRUE)
sj_long <- rbindlist(lapply(sj_files, function(f) {
  d <- fread(f, sep='\t', header=FALSE,
    col.names=c('chr','start','end','strand','motif','annot','unique','multi','overhang'))
  d$coord.intron <- paste(d$chr, d$start, d$end, sep=':')
  d$sample <- sub('_SJ.out.tab$', '', basename(f))
  d[, .(coord.intron, sample, unique)]
}))
sj <- dcast(sj_long, coord.intron ~ sample, value.var='unique', fill=0)

features <- lapply(c('SE','A5SS','A3SS','MXE','RI'), function(ev)
  read.table(sprintf('events_%s.txt', ev), header=TRUE, sep='\t'))
names(features) <- c('SE','A5SS','A3SS','MXE','RI')
seurat_obj <- readRDS('cells.rds')
pheno <- seurat_obj@meta.data; pheno$sample.id <- rownames(pheno)

marvel <- CreateMarvelObject(
  SpliceJunction=sj, SplicePheno=pheno, SpliceFeature=features,
  IntronCounts=read.table('intron_counts.tsv', header=TRUE, sep='\t'),
  GeneFeature=read.table('gene_features.tsv', header=TRUE, sep='\t'),
  Exp=read.table('tpm.tsv', header=TRUE, sep='\t'),
  GTF=gtf)
marvel <- CheckAlignment(marvel, level='SJ')
marvel <- ComputePSI(marvel, CoverageThreshold=10, EventType='SE')
```

`Exp` must retain its `gene_id` column: do not read it with `row.names=1`. Use the same chromosome convention in `coord.intron` and the rMATS `tran_id` values.

## Input contract and RI

RI additionally requires an `IntronCounts` matrix keyed by `coord.intron`. On MARVEL 2.0.5 it also requires at least two threads and an explicit read length; otherwise it can fail with `argument is of length zero`.

```r
marvel <- ComputePSI(marvel, CoverageThreshold=10, EventType='RI',
                     thread=2, read.length=100)
```

Set `read.length` to the actual aligned read length. A coverage threshold above the data depth produces all-NA PSI. Run one event class per call, then check alignment:

```r
marvel <- CheckAlignment(marvel, level='splicing')
marvel <- CheckAlignment(marvel, level='gene')
marvel <- TransformExpValues(marvel, offset=1, transformation='log2', threshold.lower=1)
neurons <- pheno$sample.id[pheno$cell.type == 'neuron']
glia <- pheno$sample.id[pheno$cell.type == 'glia']
marvel <- AssignModality(marvel, sample.ids=neurons, min.cells=5, seed=1)
marvel <- CompareValues(marvel, cell.group.g1=neurons, cell.group.g2=glia,
  min.cells=5, method='wilcox', method.adjust='fdr', level='splicing', event.type='SE')
res <- marvel$DE$PSI$Table[['wilcox']]
```

MARVEL stores group means and `mean.diff = mean.g2 - mean.g1` as PSI times 100. `method='dts'` requires the separately installed `twosamples` package. A table with no minus-strand event can trigger MARVEL 2.0.5 preprocessing defects; retain a representative event of both strands when preparing an event table.
