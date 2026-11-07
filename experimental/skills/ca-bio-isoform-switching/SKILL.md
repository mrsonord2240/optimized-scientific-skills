---
name: ca-bio-isoform-switching
description: Use when testing differential transcript usage and interpreting isoform switch consequences.
category: Data Analysis
license: MIT
author: GPTomics
---

# Isoform Switching

Choose one route, read it, and run its command before writing new code. Paths are relative to this Skill.

| Data or question | Read |
|---|---|
| Salmon quantifications and switch consequences | `routes/switches.md` |
| Gene and transcript stage-wise DTU | `routes/manual-dtu.md` |
| Gibbs samples and uncertainty-aware transcript expression | `routes/swish.md` |
| Functional annotation and NMD consequences | `routes/consequences.md` |
| Long-read or single-cell counts | `routes/counts.md` |
| Related analyses | `routes/related.md` |

DTU is a change in transcript proportions; swish tests expression rather than DTU. Match transcript IDs, samples and annotation; one replicate per group cannot support a statistical test.
Source checked versions: IsoformSwitchAnalyzeR 2.6.0, DEXSeq 1.52.0, satuRn 1.14.0, fishpond 2.12.0.
