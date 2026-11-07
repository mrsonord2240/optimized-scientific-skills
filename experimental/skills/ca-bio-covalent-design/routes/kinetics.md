# Kinetic interpretation

```python
import pandas as pd
data = pd.read_csv('fitted_kinetics.csv')
if data[['kinact_s_inv', 'Ki_M']].isna().any().any() or (data[['kinact_s_inv', 'Ki_M']] <= 0).any().any():
    raise ValueError('kinact and Ki must be positive and nonmissing')
data['kinact_over_Ki_M_inv_s_inv'] = data['kinact_s_inv'] / data['Ki_M']
data.to_csv('kinetic_efficiency.csv', index=False)
```

Supply positive, already-fitted kinact in s^-1 and Ki in M for an appropriate irreversible two-step model; reject zero, negative or missing values before calculating. This calculation does not fit a kinetic model. Compare matched assays and report fitted values and conditions. Reversible systems may need other models; GSH half-life depends on the compound and assay. See `references/kinetics.md` and `references/troubleshooting.md`.

Done when kinetic efficiency, assay conditions and model assumptions are documented.
