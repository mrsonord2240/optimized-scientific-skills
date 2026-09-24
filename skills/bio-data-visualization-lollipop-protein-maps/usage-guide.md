# Lollipop / Needle Protein Maps - Usage Guide

Use this Skill for a single-gene mutation map where recurrence, protein position, and domain context matter. It provides maftools for a standard MAF plot, trackViewer for a verified custom domain map, and g3viz for an HTML supplement.

Example prompts:

- "Plot TP53 from this MAF, state whether heights are mutation rows or unique samples, and label R175, R248, and R273."
- "Compare TP53 in Luminal and Basal samples, retaining the variant caller's protein-change column."
- "Use a TP53 P04637 trackViewer map, distinguish current UniProt features from the conventional 102--292 DNA-binding-core range, and warn about non-parseable HGVSp strings."
- "Create a g3viz HTML supplement and save it with htmlwidgets."

Before requesting a figure, provide the MAF, its protein-change column, and the transcript/protein convention used upstream. The agent-facing instructions, isoform constraints, runnable code, and failure checks live in `SKILL.md`.

## Related Skills

- data-visualization/oncoprint-mutation-matrices
- variant-calling/variant-annotation
- data-visualization/color-palettes
