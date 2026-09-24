# Distribution Plots - Usage Guide

## Overview

Use this Skill to replace bars of means with an honest per-group distribution view. `SKILL.md` contains the N-based choice, bandwidth guard, package compatibility, and failure modes; use its standalone R example for a current ggplot2 raincloud.

## Prerequisites

R: `ggplot2`, `ggbeeswarm`, `ggdist`, `lvplot`; add `introdataviz` only when a guarded split violin is needed. Python: `seaborn`, `matplotlib`, and optionally `ptitprince` (see its compatibility limit in `SKILL.md`).

## Example Prompts

- "Replace this bar plot with a ggdist raincloud showing each individual point and N per treatment arm."
- "Make a quasirandom beeswarm of expression across six cell types with n=15 per group."
- "Make a split violin of gene X by cluster for ordered Control and Treatment, omitting clusters without both conditions or fewer than 30 cells per condition."
- "Use the safe SJ-to-nrd0 bandwidth guard for these zero-inflated clusters."
- "Use a letter-value plot for N=2000 per group instead of a standard boxplot."

## Related Skills

- data-visualization/statistical-annotation - p-value brackets between groups
- data-visualization/color-palettes - CVD-safe categorical palettes
- data-visualization/ggplot2-fundamentals - underlying grammar
- single-cell/markers-annotation - stacked / split violin for scRNA
- clinical-biostatistics/effect-measures - effect size to report alongside distribution
