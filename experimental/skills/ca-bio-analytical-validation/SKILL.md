---
name: ca-bio-analytical-validation
category: Data Analysis
description: Use when stating, checking or designing the limit of detection of a ctDNA sequencing assay, including how many genome equivalents a target VAF needs.
tool_type: python
primary_tool: scipy
license: MIT
author: GPTomics
---

# ctDNA Analytical Validation and Detection Limits

Pick the one row that matches your question. Read that file and only that file. Run its command before writing any code of your own.
Paths are relative to this Skill's directory.

| Your question | Read |
|---------------|------|
| How many genome equivalents, or how much input, for a target VAF; chance the variant is in the tube at all | `routes/ge-poisson.md` |
| LoD95 (and LoB) from a replicated dilution series | `routes/lod95-probit.md` |
| Best-case sensitivity of a multi-variant panel with a k-of-N rule | `routes/panel-integrated-lod.md` |
| Check a claim such as "detects 0.1% VAF", or choose between LoD and LoQ | `routes/audit-claim.md` |
| Which reference material to validate against (HCC1395, SEQC2 Sample A/B) | `routes/reference-materials.md` |

Rules on every route:

- A VAF is a sensitivity claim only with its input mass in genome equivalents (GE), the unique-molecule (consensus) depth and the replicate detection rate. A bare "0.1% VAF" is undefined.
- GE convention: 3.3 pg per haploid genome, so 303.03 GE per ng. Another convention must be stated and passed as `--ge-per-ng`.
- The Poisson and k-of-N results are sampling only: they leave out recovery, consensus depth, background error, false positives and the calling rule. Never call them an assay LoD or LoD95.
- An assay LoD95 comes only from replicated full-workflow detections at three or more VAF levels that bracket 95%, reported with its confidence interval.
- VAF is a fraction everywhere: 0.1% is 0.001.
- Sequencing depth past sampling saturation re-reads the same molecules. More plasma or better conversion helps; more reads do not.
- Neighbouring Skills: ctdna-mutation-detection (low-VAF calls), longitudinal-monitoring, tumor-fraction-estimation, experimental-design/multiple-testing, clinical-biostatistics/power-and-sample-size.

Tested with numpy 2.5.3, scipy 1.18.1, statsmodels 0.15.0, pandas 3.0.6.
