# Interpretation


| Metric | Recommendation |
|---|---|
| Cells supporting an event | At least 50 for per-cell conclusions |
| Junction reads per event per cell | 5-10 for stable PSI; <=1 is unreliable |
| Per-cell PSI variance | Look for low within-cluster and large between-cluster separation, but inspect coverage first |
| Pseudobulk replicates | >=3 samples per group; >=50 cells per cluster is a practical minimum |
| Library choice | Full-length plate or long-read for transcriptome-wide AS; 3' chemistry for APA |

MARVEL modality labels are descriptive, not proof of mechanism: included (PSI near 1), excluded (near 0), bimodal (mixed states or bursting), middle (often technical mixture/low coverage), and multimodal. Confirm unusual modality with full-length data and doublet checks.

Beta-binomial or Dirichlet-multinomial models help with sparse overdispersed counts; they do not make unsupported 3' per-cell PSI reliable. Pseudobulk trades within-cluster heterogeneity for statistical power.

## Reconcile Disagreement

| Pattern | Likely cause | Next action |
|---|---|---|
| MARVEL significant, BRIE2 not | BRIE2 is more conservative about per-cell uncertainty | Use MARVEL for cell-type comparisons; retain BRIE2 for within-cluster uncertainty |
| BRIE2 significant, MARVEL not | Smooth cell-state association, not a discrete boundary | Test the trajectory with Psix |
| SpliZ significant, MARVEL not | Novel or differently represented junction structure | Inspect junctions manually |
| Sierra significant, MARVEL not | APA and splicing are different biology | Report them separately |
| Pseudobulk significant, per-cell method not | Aggregate power exceeds per-cell power | Report a cluster-level result |
