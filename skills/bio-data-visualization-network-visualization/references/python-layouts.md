# Python Static Layouts

Read this file when selecting a NetworkX layout, producing a static PPI figure, or rendering a signed
directed regulatory network. The main `SKILL.md` owns interpretation rules and thresholds.

## Reproducible Layout Comparison

The moved layout recipe is now a complete CLI with imports, GraphML input, finite-coordinate checks,
fixed seeds for stochastic methods, community color, degree size, and adaptive labels:

```bash
python scripts/layouts.py network.graphml --output-dir out/layouts
```

It renders spring/Fruchterman-Reingold, Kamada-Kawai, circular, spectral, and NetworkX's native
ForceAtlas2. Use `--layouts spring,forceatlas2` to select a subset. ForceAtlas2 is available directly
in networkx 3.4+; no `fa2_modified` package or Gephi round trip is required for this recipe.

- Spring is a general default for modest undirected graphs.
- Kamada-Kawai often separates small dense graphs cleanly.
- Circular emphasizes symmetry rather than inferred distance.
- Spectral can expose connectivity structure but does not turn coordinates into measurements.
- ForceAtlas2 is useful for larger hub-spoke networks, subject to runtime and legibility checks.

## Standard Static PPI

The complete static example owns the visible-width helper and the shared node/legend color mapping:

```bash
python examples/network_plots.py --graphml network.graphml --output-dir out/static --label TP53
```

Replace its deterministic demo graph with a real loader while keeping the rendering functions. The
example computes the adaptive top-k set once and uses it in every view; repeat `--label` for
prespecified genes of interest. The hub panel emphasizes its top five nodes without changing which
labels are eligible.

## Signed Directed Regulatory Network

The directed recipe uses Graphviz `dot` for hierarchy, arrowheads for direction, and separate colors
for `+` activation and `-` repression:

```bash
python scripts/directed_network.py grn.graphml --output out/grn.png --sign-attribute sign
```

This requires both `pydot` and the Graphviz `dot` executable. Do not silently coerce an undirected
network: the script rejects it. It also raises a `ValueError` naming unknown sign values rather than
writing a plausible node-only figure. Normalize values such as `activation` and `repression` to `+`
and `-` before rendering.

## Shared Layout for Two Conditions

Build a union graph, compute one position mapping, and reuse it:

```python
union = nx.compose(control, treatment)
pos = nx.spring_layout(union, seed=42)
```

Pass `pos` unchanged to both draws and use the same color, size, and width scales.
