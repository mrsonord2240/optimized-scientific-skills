# Interactive HTML and Cytoscape

Read this file for PyVis HTML or py4cytoscape desktop automation.

## PyVis

Run the complete example:

```bash
python examples/interactive_network.py --graphml network.graphml --output-dir out/interactive
```

The basic route calls `from_nx(graph.copy())` because PyVis rewrites edge attributes such as `weight`
into display fields. The styled route defines its palette locally and adds nodes and edges explicitly.
Both create `Network(directed=graph.is_directed())`, so a regulatory `DiGraph` keeps arrowheads.
Compute graph statistics and communities before handing any graph copy to PyVis. The save helper
removes PyVis 0.3.2's duplicated template headings and writes one escaped title.

PyVis wraps vis.js and is useful for a supplementary HTML view, not a static journal figure. Treat
about 2,000 nodes as a starting limit and inspect actual responsiveness and file size.

## Cytoscape Automation

Start Cytoscape Desktop, then run:

```bash
python examples/cytoscape_automation.py --graphml network.graphml --output-dir out/cytoscape
```

The example creates the mapped `degree`, `gene_type`, `score`, and capped `display_label` columns before
upload. It labels the same adaptive degree-ranked subset used by the static workflow, applies the
force-directed layout and final style, scales the layout 2.4-fold along the landscape export's x axis,
spaces the capped hubs evenly around the force layout's perimeter with their labels anchored outward,
then calls `fit_content()`. The style uses an
explicit 18-point label font and an 18-56 node-size range, so the capped labels remain readable instead
of being covered by oversized hubs. With py4cytoscape 1.13,
`set_node_shape_mapping` is discrete and does not accept `mapping_type`. Exports use resolved absolute
paths and `overwrite_file=True`; failures print a traceback and return a nonzero exit. For a portrait
canvas, replace the x-axis spacing with y-axis spacing before the final fit.

The shipped style is for an undirected PPI and deliberately has no target arrows. For a regulatory
network, explicitly set `EDGE_TARGET_ARROW_SHAPE` and map the sign attribute to edge colors, or use the
tested static directed recipe in `references/python-layouts.md`.
