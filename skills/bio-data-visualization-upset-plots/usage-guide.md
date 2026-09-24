# UpSet Plots - Usage Guide

Use this Skill to compare membership across four or more gene, peak, variant, or feature sets when a Venn diagram is no longer readable. It defaults to ComplexUpset; use UpSetR only to reproduce a legacy figure.

## Prerequisites

```r
install.packages("ComplexUpset")
# Legacy figures only: install.packages("UpSetR")
```

```bash
pip install "upsetplot==0.9.0" "pandas>=2.2,<3"
```

The full input checks, executable examples, limitations, and export settings are in [SKILL.md](SKILL.md).

## Example prompts

- “Make an UpSet plot of six gene sets, after trimming IDs and stopping if a set is empty. Sort by intersection size.”
- “Show only intersections shared by at least two sets, and export a Type-42 PDF.”
- “Add log2FC and significant metadata to a ComplexUpset plot.”
- “Highlight SetA-and-SetB only; create separate figures if I need several highlights.”
- “Reproduce this UpSetR figure without dropping any of my seven sets.”
- “Use upsetplot with pandas 2.2–2.x, exact-membership styling, and manual count labels.”

## Related Skills

- data-visualization/heatmaps-clustering
- pathway-analysis/go-enrichment
- differential-expression/de-results
