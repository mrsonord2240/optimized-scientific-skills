---
name: bio-data-visualization-network-visualization
description: Visualize biological networks (PPI, gene-regulatory, co-expression, pathway) with reproducible NetworkX layouts, degree and community encodings, signed directed edges, R hierarchical edge bundling, PyVis HTML, and Cytoscape automation. Use for static publication figures, interactive exploration, or Cytoscape export.
tool_type: python
primary_tool: NetworkX
---

# Network Visualization

Select a layout that matches the network and question, encode node and edge attributes explicitly,
and treat force-directed positions as drawing artifacts rather than biological coordinates. Use
NetworkX plus matplotlib for static figures, PyVis for modest interactive networks, ggraph for R
workflows and hierarchical edge bundling, or Cytoscape for desktop compositing.

## Version Compatibility and Prerequisites

The executable recipes were checked with networkx 3.6.1, matplotlib 3.11.2, pyvis 0.3.2,
py4cytoscape 1.13.0, Cytoscape 3.10.4, pydot 4.0.1, Graphviz 13.1.2, igraph 2.3.0,
ggraph 2.2.2, and ggplot2 4.0.3. NetworkX's native ForceAtlas2 layout requires networkx 3.4+.

```bash
pip install "networkx>=3.4" matplotlib numpy pandas pyvis py4cytoscape pydot
```

Install the Graphviz system executable as well as `pydot` before using the directed `dot` recipe.
Cytoscape Desktop must be running before the py4cytoscape example starts.

```r
install.packages(c('igraph', 'ggraph', 'ggplot2'))
```

If a signature differs, inspect it instead of retrying unchanged: `help(module.function)` in Python,
or `packageVersion('<pkg>')` and `?function_name` in R.

## The Central Constraint: Layout Is Not Biology

A force-directed layout minimizes a drawing objective involving attraction and repulsion. A node's
position therefore depends on algorithm, initialization, iteration count, and parameters; it is not
a measured biological coordinate. Set a seed for every stochastic layout and interpret only explicit
graph properties such as edge existence, degree, direction, sign, or an independently computed module.
Do not infer that one module is biologically closer to another because the drawing places them nearby.

## Decision Tree

| Network or question | Route | Why |
| --- | --- | --- |
| Generic PPI below about 500 nodes | Fruchterman-Reingold or Kamada-Kawai; read `references/python-layouts.md` | General-purpose separation |
| Scale-free PPI above about 500 nodes | Native NetworkX ForceAtlas2; read `references/python-layouts.md` | Hub-aware force-directed layout |
| Gene-regulatory network | Graphviz `dot` plus signed arrow styling; read `references/python-layouts.md` | Direction and regulatory sign remain visible |
| Pathway or manually curated signaling map | Cytoscape; read `references/interactive-and-cytoscape.md` | Desktop layout and compositing |
| Hierarchical many-to-many relations | ggraph connection bundling; read `references/r-layouts.md` | Routes relations along a known hierarchy |
| Small dense graph | Circular layout; read `references/python-layouts.md` | Preserves symmetry without implying distance |
| Connectivity-only question | Adjacency matrix heatmap | Avoids layout artifacts entirely |
| Two conditions | Compute one union layout and reuse it | Prevents layout changes from masquerading as biology |
| Too large for a legible node-link view | Filter or aggregate with a stated rule, or use an adjacency summary | A raw hairball is not informative |

## Reference Files

- Read `references/python-layouts.md` for executable static layouts, adaptive visual mappings, signed
  directed networks, and the moved NetworkX recipes.
- Read `references/r-layouts.md` only for igraph/ggraph layouts or hierarchical edge bundling.
- Read `references/interactive-and-cytoscape.md` for PyVis behavior, directed HTML, and Cytoscape
  styling/export automation.

## Core Workflow

1. Load an edge list, GraphML, SIF, or interaction-database result; retain direction, sign, and weight.
2. Report node count, edge count, directedness, connected components, and missing attributes.
3. Choose a route from the decision tree and record the algorithm, seed, and relevant parameters.
4. Compute degree or centrality and communities independently of the renderer.
5. Encode node size, node color, edge width, direction, and sign from named attributes.
6. Label a ranked subset rather than applying a fixed absolute degree cutoff.
7. Export to an explicit path and inspect the rendered output, not merely the process exit code.
8. State that positions from stochastic force-directed layouts are not biological measurements.

For the standard PPI demonstration:

```bash
python examples/network_plots.py --graphml network.graphml --output-dir out/static
```

This example uses one discrete color mapping for both nodes and legend swatches, ranks hubs for
labels in every view, treats missing edge weights as 1, and rescales widths to 0.5-4 points. Repeat
`--label GENE` to retain prespecified genes of interest in addition to the adaptive ranked set.

## Comparing Conditions

Construct the union graph, compute `pos` once, and pass that same mapping to both renders. Keep visual
scales identical. A node absent from one condition may still retain its union coordinate when shown as
missing or muted. Never run `spring_layout` separately and interpret movement as biological change.

## Per-Method Failure Modes

### Force-directed positions interpreted as biology

**Trigger:** A conclusion uses visual distance or an apparent path through empty space.

**Mechanism:** Positions minimize a drawing objective and change with initialization.

**Fix:** Base conclusions on edges, graph statistics, or orthogonal biological evidence. Use the same
seed only for reproducibility, not as evidence that a position is meaningful.

### Layout differs across runs

**Trigger:** A stochastic layout is called without a seed.

**Symptom:** Re-running the same graph produces a visibly different figure.

**Fix:** Use `seed=42` for NetworkX spring/ForceAtlas2 or `set.seed(42)` before R layouts, and report it.

### Two conditions use independent layouts

**Trigger:** Each condition calls a layout function separately.

**Symptom:** Stable nodes appear to move.

**Fix:** Compute positions on the union graph and reuse them.

### Labels are absent or overwhelm the graph

**Trigger:** A fixed rule such as `degree >= 10` is applied across graphs of different density.

**Symptom:** A sparse real network has no labels while a dense synthetic graph has dozens.

**Fix:** Rank by degree or the relevant centrality and label a capped top-k set (threshold below), plus
pre-specified genes of interest.

### Edge widths are invisible or crash on unweighted input

**Trigger:** Direct `G[u][v]['weight']` access or raw sub-point confidence values.

**Symptom:** `KeyError: 'weight'`, or every edge looks equally thin.

**Fix:** Read `data.get('weight', 1.0)` and linearly rescale the observed range to 0.5-4 points. Use a
single mid-range width when all values are identical.

### PyVis changes graph attributes

**Trigger:** `net.from_nx(G)` receives the caller's graph.

**Mechanism:** PyVis converts attributes such as `weight` into display fields in place.

**Fix:** Pass `G.copy()`, or add nodes and edges explicitly. Compute communities before conversion.

### Direction or regulatory sign disappears

**Trigger:** A directed graph is rendered by an undirected PyVis network or by Cytoscape's default
style, whose target arrow shape is `NONE`.

**Fix:** Use `Network(directed=True)` for PyVis. For static GRNs, use the signed Graphviz recipe. If
using Cytoscape, explicitly map the target arrow shape and sign colors. Validate the sign domain
before drawing so an unknown value cannot silently disappear.

### Cytoscape silently keeps default styling

**Trigger:** A mapping refers to a missing node column, or a broad exception hides an API error.

**Fix:** Populate mapped attributes before upload, allow the traceback to surface, and read back at
least one visual property. The shipped example creates `degree` and uses the py4cytoscape 1.13 shape
signature without a nonexistent `mapping_type` argument.

### Dense hairball

**Trigger:** Thousands of edges are drawn without filtering or a hierarchy.

**Fix:** Apply a declared confidence filter, use hierarchical bundling only when a real hierarchy
exists, or present an adjacency summary. Bundling is not valid without a backbone hierarchy.

## Quantitative Thresholds

| Decision | Operational value | Interpretation |
| --- | --- | --- |
| Hub labels | top `min(N, max(15, ceil(0.02N)))`, capped at 30 | Adaptive across sparse and dense graphs |
| Static edge width | 0.5-4 points after normalization | Keeps weak and strong edges visible |
| PyVis HTML | roughly below 2,000 nodes | Inspect browser responsiveness and file size |
| Consider hierarchical bundling | above about 5,000 edges, only with a hierarchy | Reduces clutter without inventing structure |
| Random seed | always record one for stochastic layouts | Reproducibility, not biological meaning |

Treat the numeric size limits as practical starting points, not universal scientific cutoffs.

## Common Errors

| Error or symptom | Cause | Solution |
| --- | --- | --- |
| `NameError: np is not defined` | Fragment copied without imports | Run the shipped script rather than a detached fragment |
| `Node ... has no position` | Positions computed on only one condition | Compute positions on the union graph |
| `KeyError: 'weight'` | Unweighted input | Use a default and normalize |
| PyVis edge weights disappear | `from_nx` mutated the graph | Pass a copy |
| GRN has no arrows | Renderer created as undirected or Cytoscape default | Enable direction explicitly |
| Cytoscape export already exists | Overwrite not enabled | Use `overwrite_file=True` and an absolute path |
| Graphviz layout fails | `pydot` or the `dot` executable is absent | Install both prerequisites |
| HTML stalls | Interactive graph exceeds browser capacity | Filter, aggregate, or use a non-node-link summary |

## References

- Csardi G, Nepusz T. 2006. The igraph software package for complex network research.
- Fruchterman TMJ, Reingold EM. 1991. Graph drawing by force-directed placement.
- Hagberg A, Schult D, Swart P. 2008. Exploring network structure, dynamics, and function using NetworkX.
- Holten D. 2006. Hierarchical edge bundles: visualization of adjacency relations in hierarchical data.
- Jacomy M, et al. 2014. ForceAtlas2, a continuous graph layout algorithm for handy network visualization.
- Shannon P, et al. 2003. Cytoscape: a software environment for integrated models of biomolecular interaction networks.

## Related Skills

- gene-regulatory-networks/coexpression-networks - Build the network to visualize
- database-access/interaction-databases - Fetch PPI data
- data-visualization/multipanel-figures - Combine network with other plots
- data-visualization/color-palettes - Choose accessible community colors
- single-cell/cell-communication - Visualize cell-cell interaction networks
