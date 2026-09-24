# Tree Visualization - Usage Guide

## Overview

This skill draws, styles, annotates, and exports phylogenetic tree figures, and chooses the drawing tool. A tree figure is an argument, not a neutral picture: layout geometry, node ordering, root placement, and which support value is shown are all supplied by whoever draws the tree. The skill makes those choices conscious and declares them in the caption. The agent's method, thresholds, and failure modes are in `SKILL.md`; this guide is for choosing the skill.

## Quick Start

Tell your AI agent what you want to do:
- "Draw this Newick tree as a phylogram and save it as a vector PDF"
- "Show bootstrap support at the internal nodes and label it as bootstrap in the caption"
- "Color the branches for one clade red and ladderize for legibility"
- "My tree has 800 tips and the labels collide - what layout should I use?"
- "I have a BEAST MCC tree with HPD bars - how do I draw the node-age uncertainty?"

## Example Prompts

### Drawing and exporting
> "Read tree.nwk, ladderize it, draw it as a phylogram with a title naming the branch-length unit, and export to tree.pdf as vector"
> "Give me an ASCII diagram of this tree so I can sanity-check the topology in a log"

### Layout choice for size
> "This tree has ~600 tips and rectangular labels are an unreadable band - recommend and produce a circular layout with radial labels"
> "I want metadata rings (host species, sampling year) around a circular tree of 400 genomes"

### Support and honesty
> "Show the SH-aLRT/UFBoot dual support at each node and state which is which in the legend"
> "Collapse every node below 70% bootstrap into a polytomy and say so in the caption"

### Annotated Bayesian trees
> "Draw this BEAST chronogram with the 95% HPD bars on node ages and the posterior probabilities labeled"
> "Bio.Phylo dropped my HPD intervals - route this through treeio and ggtree instead"

## Related Skills

- tree-io - parsing and preserving BEAST/MrBayes/IQ-TREE annotations so they survive into the figure
- tree-manipulation - rooting, ladderizing, pruning, and collapsing low-support nodes before drawing
- data-visualization/ggplot2-fundamentals - the grammar, themes, and ggsave vector export underlying ggtree
- data-visualization/multipanel-figures - composing a tree plus aligned metadata panels into one figure
