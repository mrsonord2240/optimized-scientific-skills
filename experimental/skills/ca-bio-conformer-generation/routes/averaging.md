# Ensemble-weighted property

```python
import pandas as pd
from scripts.gen_conformers import boltzmann_average
data = pd.read_csv('conformer_properties.csv')
value = boltzmann_average(data['value'].to_numpy(), data['energy_kcal_mol'].to_numpy(), T=300.0)
pd.DataFrame([{'boltzmann_average': value, 'temperature_K': 300.0}]).to_csv('average.csv', index=False)
```

The helper is `scripts/gen_conformers.py`. Supply one matched value and energy per conformer in kcal/mol; change temperature only as justified. See `references/averaging.md`. Raw gas-phase force-field minima are exploratory surrogates and do not establish thermodynamic populations.

Done when the property summary and population assumptions are reported.
