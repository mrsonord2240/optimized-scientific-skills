# Network Visualization - Usage Guide

## Overview

Use this Skill to turn protein-interaction, gene-regulatory, co-expression, or pathway graphs into
reproducible static figures, interactive HTML, or Cytoscape exports. The Skill selects a layout,
encodes graph attributes, protects directed and weighted data, and prevents visual layout from being
misread as a biological coordinate. Installation, thresholds, failure modes, and executable routes
live in `SKILL.md`; method details live in its indexed reference files.

## Example Prompts

### Standard PPI

> "Render this STRING PPI with a fixed spring-layout seed. Size nodes by degree, color them by
> community, normalize confidence to visible edge widths, and label a capped top-k set of hubs."

### Layout Comparison

> "Compare spring, Kamada-Kawai, circular, spectral, and native NetworkX ForceAtlas2 layouts for this
> GraphML network. Keep stochastic seeds fixed and do not interpret visual distance as biology."

### Directed Regulatory Network

> "Render this directed GRN with a Graphviz hierarchy, arrowheads, and distinct colors for activation
> and repression from the `sign` edge attribute."

### Comparing Conditions

> "Render control and treatment networks side by side. Compute positions once on the union graph and
> reuse the exact layout and visual scales in both panels."

### Interactive Supplement

> "Create a directed PyVis HTML supplement without mutating my NetworkX edge attributes; include
> degree-sized nodes, module colors, tooltips, and physics controls."

### Cytoscape Pipeline

> "Send this attributed NetworkX PPI to Cytoscape, apply degree and confidence mappings, and export
> overwriteable PDF and PNG files to an explicit output directory."

### Hierarchical Relations

> "Use ggraph connection bundling for these relations along the supplied biological hierarchy, after
> verifying that every relation endpoint exists in the hierarchy."

## Related Skills

- gene-regulatory-networks/coexpression-networks - Build the network
- database-access/interaction-databases - Fetch PPI data
- data-visualization/multipanel-figures - Combine with other plots
- data-visualization/color-palettes - Choose an accessible community palette
- single-cell/cell-communication - Visualize cell-cell interactions
