# Color Palettes for Scientific Visualization - Usage Guide

## Overview

Use this Skill to choose sequential, diverging, cyclic, or categorical color encodings and to verify grayscale and color-vision-deficiency behavior. The decision table, thresholds, executable checks, and failure modes live in `SKILL.md`.

## Prerequisites

```r
install.packages(c('viridis', 'RColorBrewer', 'scico', 'ggsci', 'colorspace'))
```

```bash
pip install matplotlib cmcrameri colorspacious scipy
```

## Example Prompts

- "Pick a sequential perceptually uniform colormap for an expression heatmap and verify its luminance remains monotonic in grayscale."
- "Use a diverging Crameri vik palette for a log-fold-change heatmap with symmetric bounds at +/- the 99th percentile of |LFC|; anchor zero at the palette's near-neutral centre."
- "Assign one stable, named Okabe-Ito color to each of these seven cell types, reserve light grey for ambient/unassigned, and add a redundant marker if CVD simulation shows a collision."
- "Find every plot in this notebook that uses `cmap='jet'` or rainbow and replace sequential uses with viridis or batlow and diverging uses with vik or roma."
- "Simulate the current palette with `colorspace::deutan()` and `protan()`, report the minimum pairwise distance, and inspect whether the categories remain distinguishable."

For implementation details, see `SKILL.md` sections **Palette Type by Data Type**, **CVD Simulation**, **Grayscale Monotonicity Test**, and **Diverging Palette Setup**.

## Related Skills

- data-visualization/heatmaps-clustering - Robust diverging bounds for heatmaps
- data-visualization/volcano-and-ma-plots - Okabe-Ito Up/Down/NS conventions
- data-visualization/ggplot2-fundamentals - Applying palettes in ggplot2 scales
- data-visualization/dimensionality-reduction-plots - Categorical palettes for cluster labels
