Moved out of `SKILL.md`: read when running PoPS on MAGMA output.

## PoPS Polygenic Priority Score

**Goal:** Add a distance-orthogonal similarity-based prior to ranked gene candidates.

**Approach:** Run MAGMA genome-wide to produce gene Z; feed gene Z plus a curated gene-feature matrix (pathway membership + co-expression + PPI) to PoPS ridge (L2-penalized) regression. Per-gene priority scores are produced; per-locus relative ranking is informative.

```bash
# PoPS requires the gene-feature matrix and MAGMA gene Z output
# Download features and gene_annot from FinucaneLab/pops releases

python pops.py \
    --gene_annot_path gene_annot.txt \
    --feature_mat_prefix PoPS_features_full \
    --control_features_path control.features \
    --num_feature_chunks 10 \
    --magma_prefix gene_step \
    --out_prefix pops_out

# Output pops_out.preds: per-gene priority score
# Output pops_out.coefs: per-feature ridge coefficients (interpretation)
```

PoPS is biology-agnostic; the feature matrix encodes biology. Bias in the features (e.g. cancer-pathway-heavy gene sets for a non-cancer trait) propagates to the output; verify feature coverage matches the trait.

**Windows: rename MAGMA's output before running PoPS.** MAGMA 1.10 on Windows writes `<magma_prefix>.genes.out.txt` (extra `.txt`), but `pops.py` hard-codes `<magma_prefix>.genes.out` and raises `FileNotFoundError` if the exact name is missing -- following this section verbatim on Windows fails with no explanation (verified end-to-end: PoPS confirmed working immediately after `cp <prefix>.genes.out.txt <prefix>.genes.out`; corroborated independently in this Skill's tooling environment, `mendelian-randomization-analyst/TOOLS.md`). `examples/pops_run.py` wraps this exact PoPS invocation and performs the rename automatically (no-op on Linux/Mac, where MAGMA writes `.genes.out` directly) -- prefer it over calling `pops.py` directly on Windows.

**PoPS needs genome-wide, multi-chromosome MAGMA input to be meaningful.** With only one chromosome (or a single locus) of genes, PoPS's held-out-chromosome ridge CV has no held-out fold to validate against, `SELECTED_CV_ALPHA` saturates at its maximum, and every `PoPS_Score` collapses toward 0 -- the code path still runs to completion (exit 0, real output files), but the scores carry no signal at that scale (verified on a real 3-gene single-locus run). Do not present a locus-scale PoPS run's scores as a real result; use it only to confirm the pipeline is wired correctly, and require genome-wide MAGMA output for an actual PoPS-based effector-gene call.
