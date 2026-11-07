---
name: ca-bio-differential-expression-deseq2-basics
category: Data Analysis
description: Use when testing which genes differ between conditions from RNA-seq counts with DESeq2.
tool_type: r
primary_tool: DESeq2
license: MIT
author: GPTomics
---

# DESeq2 Differential Expression

Pick the one row that matches your data. Read that file and only that file. Run its command before writing any code of your own.
Paths are relative to this Skill's directory.

| Your data | Read |
|-----------|------|
| Single-cell counts, cells from several donors | `routes/pseudobulk.md` |
| Gene x sample count matrix, two groups | `routes/two-group.md` |
| Same, plus a batch, or the same subjects in both groups | `routes/paired-or-batch.md` |
| One factor with three or more levels | `routes/multi-level.md` |
| Two factors and a question about how they combine | `routes/interaction.md` |
| Salmon, kallisto or RSEM output, not counts | `routes/tximport.md` |
| Python only, no R | `routes/python.md` |

Rules on every route:

- Input is raw integer counts. Never normalized, log, TPM or batch-corrected values.
- Name the reference level. The scripts refuse to run without `ref=`.
- Take significance from `padj` in the results file the script writes. Look genes up in that file by name, never by row number.
- One run of one engine is the answer. Do not repeat it in edgeR or limma to compare.
- When the request sets its own cutoff or fold-change definition, apply it to the script's output. Do not swap the test.

After the run, only if the request needs it:

| Need | Read |
|------|------|
| A ranked list for GSEA, or fold changes for a plot | `routes/shrinkage.md` |
| A gene of interest shows `padj = NA` | `routes/padj-na.md` |
| Values for PCA, a heatmap or machine learning | `routes/transforms.md` |
| What to state in the write-up | `routes/report.md` |
| An error, or a result that looks wrong | `routes/errors.md` |

Tested with DESeq2 1.46.0, apeglm 1.28.0, ashr 2.2.63, PyDESeq2 0.5.4 (R 4.4).
