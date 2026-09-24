# Interpreting BAGEL2 calls (essential, neutral, tumor suppressor)

> Moved out of SKILL.md. "Below" and "Failure Modes" in the text refer to sections of SKILL.md.

## Interpret BAGEL2 Results

**Goal:** Stratify genes into essential, non-essential, and tumor-suppressor categories.

**Approach:** Apply BF threshold to classify essentials; only classify negative-BF genes as tumor suppressors when the screen design actually expects enrichment (see Failure Modes below) (verified below).

```bash
python scripts/interpret_bagel.py bayes_factor.txt --screen-type dropout -o calls.tsv
# or: from interpret_bagel import interpret_bagel   (run from scripts/)
```

`scripts/interpret_bagel.py` drops the assay-control pseudo-genes (`LacZ`, `luciferase`, `EGFP`; override with `--controls`), calls `essential` at BF > 6 (`--bf-essential`), and calls `tumor_suppressor` at BF < -6 only when `--screen-type enrichment|both`. It warns when more than 5% of genes get a tumor-suppressor call.

Verified on real HAP1 TKOv3 data (a T0-vs-T18 dropout screen): skipping the guard flags **86.5% of the genome** as tumor-suppressor, with three assay-control pseudo-genes (`LacZ`, `luciferase`, `EGFP`) as the top hits. With the guard, the default
(`screen_type='dropout'`) returns 0 tumor-suppressor calls and excludes the 3 control
pseudo-genes; explicit `screen_type='enrichment'` still surfaces the true tumor
suppressors TSC2 (most negative BF, -77.6) and TSC1 (seventh most negative, -64.9; DEPDC5, another mTOR-pathway gene, is second) (now that the control genes
that previously masked them are excluded) but raises the 86.6%-flagged warning so the
caller doesn't report it uncritically.

**Tumor suppressor identification:** Genes with significantly negative BF (e.g., <-6) in a screen actually designed to detect enrichment (drug-resistance, GoF) indicate fitness advantage from their loss, which is biologically distinct from "non-essential". **Do not call tumor suppressors from a pure dropout screen** -- see Failure Modes.
