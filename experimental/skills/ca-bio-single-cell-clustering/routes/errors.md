# Clustering problems

| Symptom | Cause | Fix |
|---------|-------|-----|
| One giant blob | Too few HVGs or `n_pcs` too low; or one real population | Raise HVGs (about 2000) and `n_pcs`, check the elbow; if markers stay uniform across the sweep it may be one population |
| Far too many clusters | Resolution too high; `n_pcs` adding noise | Lower resolution; check stability across nearby `n_pcs` (`routes/validate.md`) |
| Adjacent clusters share all top markers | Possible over-clustering | Covariate alignment, bootstrap stability, a formal split test (`routes/significance.md`); marker overlap alone cannot decide |
| Labels change between runs or versions | Leiden backend or `n_iterations` unpinned | The scripts pin `flavor='igraph', n_iterations=2, directed=False` and a seed |
| A cluster maps to one sample or lane | Batch effect | Integrate first (batch-integration Skill); inspect QC covariates |
| A stable cluster of stress or cycle genes | Technical split | Regress out or score and exclude cycle and stress genes; require non-stress markers |
| A cluster expresses two lineages' markers | Doublets | Run doublet detection before clustering |
| Marker p-values quoted as proof | Double-dipping | `routes/significance.md` |

Done when the symptom's fix has been applied and the sweep rerun.
