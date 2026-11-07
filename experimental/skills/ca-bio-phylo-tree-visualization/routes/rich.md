# Draw annotation-rich trees
~~~bash
Rscript scripts/unrooted.R tree.nwk unrooted.pdf
~~~
Done when the new unrooted vector figure is checked and branch-length units declared. An unrooted plot cannot justify basal or time-direction claims. For BEAST HPD bars, metadata rings, circular layouts, and dual support, adapt [treeio and ggtree recipes](../references/ggtree-ape-recipes.md) to the actual columns and installed versions. BEAST HPD uses read.beast and geom_range with center = height; verify one interval. The source composite example uses synthetic metadata and is not a finished real-metadata pipeline.
