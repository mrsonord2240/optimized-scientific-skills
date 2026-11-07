# Confounders, complex designs and follow-up

```python
import patsy, numpy as np
design = patsy.dmatrix('group + C(batch) + RIN', meta)
assert np.linalg.matrix_rank(np.asarray(design)) == design.shape[1]
```

Supply sample metadata as `meta`. A rank-deficient design cannot identify the requested effects. For covariate regression, multigroup approaches, thresholds, cross-tool comparison and event prioritization use `references/interpretation.md`; see `references/design.md` for design decisions.

Done when the design is identifiable and requested estimates or priorities are reported.
