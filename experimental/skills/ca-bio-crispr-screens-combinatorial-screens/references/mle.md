## Run Combinatorial Screen Analysis (MAGeCK MLE with Interaction Indicator)

**Goal:** Use MAGeCK MLE to estimate the effect of each gene independently and the additional effect when both genes are simultaneously perturbed.

**Approach:** Design matrix encodes single-A, single-B, double-AB conditions; the `interaction` column is set to 1 only for double-KO samples. The resulting beta for that column captures the extra effect beyond the sum of single-gene betas. Note: MAGeCK MLE does not natively perform a formal interaction-significance test, but the `interaction|beta` and `|fdr` columns serve as the GI estimate; for formal interaction testing, compute GI = observed_double_lfc - (single_A_lfc + single_B_lfc) explicitly (see Genetic Interaction (GI) Scoring above).

```bash
# Design matrix encoding double-KO as a separate "interaction" indicator
# Conditions: NT (control), A_KO, B_KO, A_B_KO
cat > combo_design.txt <<EOF
Samples         baseline    geneA       geneB       interaction
NT_r1           1           0           0           0
NT_r2           1           0           0           0
A_r1            1           1           0           0
A_r2            1           1           0           0
B_r1            1           0           1           0
B_r2            1           0           1           0
AB_r1           1           1           1           1
AB_r2           1           1           1           1
EOF

mageck mle \
    --count-table combo_counts.txt \
    --design-matrix combo_design.txt \
    --output-prefix combo_mle

# Output: per-gene beta scores per design column
# The "interaction" column beta captures additional joint effect beyond additive
```

**Interpretation of MAGeCK MLE output:**

| Column | Meaning |
|--------|---------|
| `geneA|beta` | Single-A effect |
| `geneB|beta` | Single-B effect |
| `interaction|beta` | Additional effect under joint perturbation beyond sum of singles |
| `interaction|p-value`, `|fdr` | Significance vs zero |

A significantly negative `interaction|beta` is synthetic lethal; positive is synthetic rescue. For formal GI hypothesis testing, prefer the explicit GI scoring approach (Genetic Interaction (GI) Scoring, above) over MAGeCK MLE interpretation, since MAGeCK MLE does not validate the additive null.

**`interaction|fdr` is not reproducible across reruns.** MAGeCK MLE's FDR is permutation-based and exposes no `--seed` flag: two consecutive `mageck mle` runs on byte-identical count/design-matrix inputs produced different `interaction|fdr` for 39 of 40 genes, while `interaction|beta` point estimates stayed bit-identical (checked on MAGeCK 0.5.9.5). Treat `interaction|beta` sign/magnitude as the stable signal; do not report a specific `interaction|fdr` value as reproducible, and use the deterministic GI z-score method above for any significance claim that must hold across reruns.
