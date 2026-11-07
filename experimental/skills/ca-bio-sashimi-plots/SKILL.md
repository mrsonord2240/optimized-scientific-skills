---
name: ca-bio-sashimi-plots
category: Data Analysis
description: Use when visualizing RNA-seq coverage and splice junctions as sashimi plots.
license: MIT
author: GPTomics
---

# Sashimi plots

Experimental. Paths are relative to this skill. Read the matching route and run its command before writing code.

| Input or goal | Route |
|---|---|
| Grouped BAM overlay or batch rMATS events | `routes/ggsashimi.md` |
| Direct rMATS event figures | `routes/rmats.md` |
| Leafcutter cluster browser | `routes/leafviz.md` |
| MAJIQ LSV viewer | `routes/majiq.md` |
| Cross-tool plots and heatmaps | `routes/jutils.md` |
| Multi-track figure | `routes/tracks.md` |

Confirm reference and contig names, verify input BAMs and rendered output, and inspect visual alignment. ggsashimi requires ggplot2 < 3.5 to avoid panel misalignment.

Tested with ggsashimi 1.1.5, rmats2sashimiplot 4.0.0, leafcutter 0.2.9, pyGenomeTracks 3.9 and Jutils 1.5. MAJIQ/VOILA commands were not run.
