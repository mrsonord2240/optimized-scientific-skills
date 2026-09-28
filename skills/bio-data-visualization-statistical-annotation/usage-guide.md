# Statistical Annotation - Usage Guide

## Overview

Use this Skill when a distribution plot needs statistical brackets and the agent must choose a design-appropriate test, adjust a declared comparison family, preserve negative results, and report effect size. The executable R and Python workflows, thresholds, correction caveats, and failure modes live in `SKILL.md`.

## Example Prompts

- "Add an adjusted Wilcoxon bracket between Control and Treatment; show the numeric adjusted p-value."
- "Run all six pairwise comparisons across four conditions with Holm adjustment and save the results table."
- "Use paired Wilcoxon for pre/post measurements keyed by `subject_id`, show the subject trajectories, and reject incomplete pairs."
- "Run Kruskal-Wallis followed by Dunn post-hoc comparisons and annotate every adjusted result."
- "These are cells nested within patients. Put the mixed-model/emmeans contrast on the violin plot rather than a cell-level test."
- "Give me the Welch and Tukey versions for an approximately normal analysis."

## Before You Ask

Provide the response column, grouping column, intended comparison family, biological replicate or subject identifier when applicable, pairing/nesting structure, and whether numeric adjusted p-values or stars are preferred. See `SKILL.md` sections "Choose the Test from the Design" and "Adjust the Intended Comparison Family" for the agent's decision flow.

## Related Skills

- data-visualization/distribution-plots - underlying distribution plot
- clinical-biostatistics/categorical-tests - categorical outcomes
- clinical-biostatistics/effect-measures - companion effect size
- experimental-design/multiple-testing - FWER and FDR choices
