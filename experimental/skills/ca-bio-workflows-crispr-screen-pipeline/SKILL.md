---
name: ca-bio-workflows-crispr-screen-pipeline
description: Use when analyzing a pooled CRISPR knockout, CRISPRi or CRISPRa screen from FASTQ or guide counts to hit genes.
category: Data Analysis
tool_type: mixed
primary_tool: MAGeCK
license: MIT
author: GPTomics
---

# CRISPR Screen Pipeline

Pick the row that matches what you have. Read that route file first (and only the earlier-stage route the order below requires), then run its command before writing any code of your own. Paths are relative to this Skill's directory.

| You have | Read |
|----------|------|
| FASTQ files and a guide library CSV | `routes/count.md` |
| A guide count table, before any hit calling | `routes/qc.md` |
| A cancer cell line screen, copy-number artifacts to remove before hit calling | `routes/cn-correction.md` |
| Two conditions (endpoint vs baseline, dropout or enrichment) | `routes/rra.md` |
| Essential-gene calls against reference gene sets | `routes/bagel2.md` |
| Time course, several conditions or batches, or a design matrix | `routes/mle.md` |
| Drug vs vehicle | `routes/drugz.md` |
| Several cell lines in one library, with a replicate map and a guide-to-gene map | `routes/jacks.md` |
| A DepMap-style panel needing copy-number-aware scores (plasmid pool, days, replicate map) | `routes/chronos.md` |
| Hit lists from two or more methods to merge | `routes/consensus.md` |
| Single-cell, combinatorial, base-editor, prime-editor or in vivo screen | `routes/branches.md` |
| A tool is missing | `references/install.md` |

Order on every screen: count, QC, copy-number correct (cancer lines), hit call, consensus. Never call hits first.

Rules on every route:

- Settle four commitments before any code, never inferred from column names alone: baseline (plasmid pool, Day-0, or vehicle), library control classes (non-targeting guides plus reference essentials and non-essentials), screen type (dropout, enrichment, FACS, drug-modifier), and copy-number profile for cancer lines. When the request or data already states one, restate it as an assumption in your reply and run the route command. Ask only for what neither states. A wrong one rescales every hit and throws no error.
- Drug screens use the vehicle arm as control, never Day-0.
- Input to every caller is the raw guide count table, never a normalized or batch-corrected one. Cancer lines: the CRISPRcleanR-corrected table.
- If known core essentials do not deplete (CEGv2 PR-AUC below 0.7), the screen failed and no novel hit is trustworthy whatever its p-value.
- Take hits from the result file a script or caller writes, looked up by gene name.

Tested with MAGeCK 0.5.9.5, BAGEL2 build 115, drugZ (2026-09 HEAD), JACKS 0.2, crispr-chronos 2.3.15 (import `chronos`), CRISPRcleanR 3.0.1 (R 4.4.3), pandas 3.0.
