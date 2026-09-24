# Computational footprint (reference for heritability-partitioning SKILL.md)

## Computational Footprint

| Method | Per-trait runtime | Hardware |
|--------|------------------|----------|
| LDSC h2 (univariate) | minutes | laptop |
| Stratified LDSC + baseline-LD | 10-20 min | laptop |
| LDSC cell-type prioritization (~200 tissues) | hours, single-threaded | server |
| HESS genome-wide local h2 | hours per chromosome | cluster |
| BOLT-REML at N = 500k | days | cluster |
| HDL.rg (genetic correlation) | seconds to minutes | laptop |
| LDAK SumHer | tens of minutes | laptop |

graphREML on biobank-scale (N > 200k) typically beats S-LDSC per-trait runtime when an LDGM panel is available. Build runtime escalates with annotation count; cell-type prioritization with ~200 tissues is the most expensive per-trait step.
