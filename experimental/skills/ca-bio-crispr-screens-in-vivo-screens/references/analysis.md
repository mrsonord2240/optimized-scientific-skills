# Analysis context and limitations
The source emphasizes that effective representation can decrease substantially through bottlenecks and that replicate variation and context-specific essentiality can complicate interpretation. In its examples, sample-level QC, matched in-vitro comparisons, per-sample MAGeCK analyses and review of meta-significant but inconsistent genes are suggested. Treat endpoint coverage thresholds and reference-set precision-recall values as context-dependent rather than universal guarantees.

## Hit Calling for In Vivo

**Goal:** Identify per-gene fitness effects despite high inter-animal variability.

**Approach:** Each animal is a "replicate" with high variance due to clonal dynamics. Use MAGeCK MLE with animal-as-batch covariate, or run MAGeCK RRA per animal and meta-analyze.

```bash
# Option A: MAGeCK MLE with batch covariate
cat > in_vivo_design.txt <<EOF
Samples         baseline    tumor   animal_2  animal_3  animal_4  animal_5
Plasmid         1           0       0         0         0         0
Animal1         1           1       0         0         0         0
Animal2         1           1       1         0         0         0
Animal3         1           1       0         1         0         0
Animal4         1           1       0         0         1         0
Animal5         1           1       0         0         0         1
EOF

mageck mle \
    --count-table in_vivo_screen.count.txt \
    --design-matrix in_vivo_design.txt \
    --output-prefix in_vivo_mle
```

**MLE determinism note (checked on MAGeCK 0.5.9.5):** `mageck mle`'s permutation-derived FDR-like p-values drift run-to-run on identical input — beta effect sizes and Wald p-values stay byte-identical, but the permutation-based FDR column does not, because `mageck mle` exposes no `--seed` option (confirmed via `mageck mle --help`). Report **Wald p-values** as the primary significance metric for MLE results; `--permutation-round` (default 2, suggested 10) improves the permutation FDR estimate's stability but does not make it deterministic.

**Per-animal RRA + meta-analysis:**

Run `mageck test` on each animal vs plasmid, then combine the per-animal `gene_summary.txt` files with Stouffer's Z, using `neg|p-value` clipped to keep `norm.ppf` finite and negated so positive z means stronger depletion:

```bash
python scripts/per_animal_meta_analysis.py   # edit animal_files at the top; writes in_vivo_meta_hits.tsv and in_vivo_meta_all.tsv
```

**Compound hit-calling threshold:** call a gene a hit when **meta-FDR < 0.05 AND >=50% of animals individually reach nominal per-animal p < 0.05** (same direction). Use per-animal **nominal** p, not per-animal FDR, for the consistency arm: per-animal FDR-correction demands more power than typical in vivo cohorts provide (5-10 animals against a gene universe of hundreds-thousands), and requiring it can silently report zero hits on a screen with an obvious, strong meta-signal.

Verified against a synthetic 6-animal / 60-gene in vivo dataset with 5 planted true hits (`meta_fdr` reaching 3e-6 via the Stouffer method above): the per-animal-**FDR**<0.05 consistency rule called **0/5** planted hits (no single animal's own MAGeCK RRA test crossed FDR<0.05 at this cohort size), while the per-animal-**nominal**-p<0.05 rule recovered **4/5** with **zero** false positives among the 55 noise genes. A gene that clears meta-FDR<0.05 but fails even the nominal-p consistency check (the 5th planted hit in that test: meta-FDR 0.0013, but nominal p<0.05 in only 1/6 animals) should be **flagged for manual review, not silently discarded** — it may reflect a real effect concentrated in a minority of animals (e.g. a tumor-specific clonal effect) or an outlier animal.


The source MLE design is saturated as written and its animal indicators cannot establish a general batch adjustment with these observations. Its Stouffer example assumes independent p-values despite a potentially shared baseline. The numerical examples and thresholds are illustrative; assess calibration for the actual dataset.
