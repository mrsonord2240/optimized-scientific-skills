# R Layouts and Hierarchical Edge Bundling

Read this file for igraph/ggraph rendering. The main `SKILL.md` owns layout interpretation and the
decision about whether a hierarchy is legitimate.

## ggraph Layouts

The repaired R layout recipe reads GraphML, fixes the random seed, and renders Fruchterman-Reingold,
Kamada-Kawai, and circular layouts as separate valid plots; it no longer ends on a dangling `+`:

```bash
Rscript scripts/ggraph_layouts.R network.graphml out/r-layouts
```

## Hierarchical Edge Bundling

Bundling requires a hierarchy plus relations between leaves. It is not a generic transform for any
edge list. The executable demonstration uses ggraph's `flare` hierarchy and its matching imports, so
every index is defined before `get_con` is called:

```bash
Rscript scripts/edge_bundling.R out/edge-bundling.png
```

For biological data, replace `flare$edges`, `flare$vertices`, and `flare$imports` together. Assert that
all relation endpoints match hierarchy vertices before plotting. The old fragment's undefined
`graph`, `from_idx`, and `to_idx` names have no place in a runnable workflow.
