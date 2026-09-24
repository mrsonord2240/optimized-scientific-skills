## ggtree + treeio Recipe (R)

A composite publication figure -- IQ-TREE support kept as columns, rooted, dual support labeled, and a metadata ring -- checked on ggtree 3.14.0 / treeio 1.30.0 / ggtreeExtra 1.16.0 / ggplot2 4.0.3:

```r
suppressPackageStartupMessages({library(treeio); library(ggtree); library(ggtreeExtra); library(ggplot2)})

iq <- read.iqtree('tree.treefile')          # keeps SH-aLRT/UFBoot as columns instead of text stuck in a label

root_keep <- function(td, outgroup) {       # treeio root() can renumber tip labels to "1","2",... -- restore by index
  r <- treeio::root(td, outgroup = outgroup, edgelabel = TRUE)
  if (all(grepl('^[0-9]+$', r@phylo$tip.label))) r@phylo$tip.label <- td@phylo$tip.label[as.integer(r@phylo$tip.label)]
  r
}
iq_r <- root_keep(iq, c('OutA', 'OutB'))    # root BEFORE labeling or coloring -- same root-first rule as the Bio.Phylo recipe in `bio-phylo-recipes.md`

p <- ggtree(iq_r, size = 0.4) +
  geom_tiplab(size = 2.6, offset = 0.003) +
  geom_nodelab(aes(label = ifelse(is.na(UFboot), '', paste0(SH_aLRT, '/', UFboot))), size = 2, hjust = 1.1, vjust = -0.5) +
  geom_treescale(width = 0.02, fontsize = 2.4)

meta <- data.frame(label = iq_r@phylo$tip.label, trait = seq_along(iq_r@phylo$tip.label))   # replace with real metadata
p2 <- p + geom_fruit(data = meta, geom = geom_tile, mapping = aes(y = label, fill = trait), pwidth = 0.06, offset = 0.08)
ggsave('fig.pdf', p2, width = 180, height = 150, units = 'mm')
```

Unrooted view with a scale bar (base R, `ape` 5.8.1 -- the fallback when ggtree's `equal_angle`/`daylight` layouts fail; an unrooted view commits to no root, so make no "basal" claim from it):

```r
library(ape)
tr <- read.tree('tree.nwk')
pdf('unrooted.pdf', width = 7, height = 7)
plot(tr, type = 'unrooted', cex = 0.6, no.margin = TRUE)
add.scale.bar(cex = 0.7, lwd = 1)               # mandatory on any phylogram; units are the file's branch-length units
dev.off()
```

Use `geom_fruit` (ggtreeExtra) for the metadata ring, not `gheatmap` -- see the Version Compatibility note in `SKILL.md` on `gheatmap`'s composite-pipeline failure. `groupOTU(tree, list(...), group_name = 'grp')` + `aes(color = grp)` colors branches by predefined clade membership rather than one MRCA at a time; checked here via `ggplot_build()` segment colours that the stem edge into a defined group is colored by that group, not left on the parent's/ungrouped default -- root first here too, for the same reason the Bio.Phylo color recipe (`bio-phylo-recipes.md`) roots first.
