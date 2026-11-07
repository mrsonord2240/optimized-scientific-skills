# R layouts and bundles

```bash
Rscript scripts/ggraph_layouts.R network.graphml out/r-layouts
```

This reads GraphML and renders seeded Fruchterman-Reingold, Kamada-Kawai and circular plots. It may overwrite outputs. The bundled hierarchy demonstration is `Rscript scripts/edge_bundling.R out/bundle.png`; it uses ggraph's example hierarchy, not user data. For another hierarchy adapt all vertex, edge and relation inputs together and validate endpoints. Bundling requires a real hierarchy and is not a generic edge-list operation. See [R details](../references/r-layouts.md).

Done when layout, graph semantics, relation mapping and rendered output have been reviewed.
