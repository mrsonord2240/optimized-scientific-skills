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
Compute graph statistics and communities before handing any graph copy to PyVis.

PyVis wraps vis.js and is useful for a supplementary HTML view, not a static journal figure. Treat
about 2,000 nodes as a starting limit and inspect actual responsiveness and file size.

## Cytoscape Automation

Start Cytoscape Desktop, then run:

```bash
python examples/cytoscape_automation.py --graphml network.graphml --output-dir out/cytoscape
```

The example creates the mapped `degree`, `gene_type`, and `score` columns before upload. With
py4cytoscape 1.13, `set_node_shape_mapping` is discrete and does not accept `mapping_type`. Exports use
resolved absolute paths and `overwrite_file=True`; failures print a traceback and return a nonzero exit.

The shipped style is for an undirected PPI and deliberately has no target arrows. For a regulatory
network, explicitly set `EDGE_TARGET_ARROW_SHAPE` and map the sign attribute to edge colors, or use the
tested static directed recipe in `references/python-layouts.md`.
