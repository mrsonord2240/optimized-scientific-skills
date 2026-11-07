# Binding and degradation measurements

```python
import pandas as pd
from scripts.cooperativity_dc50 import detect_hook, fit_dc50
data = pd.read_csv('dose_response.csv').sort_values('concentration_nM')
conc = data['concentration_nM'].to_numpy()
degradation = data['degradation_percent'].to_numpy()
hook = detect_hook(conc, degradation)
fit = fit_dc50(conc, degradation, peak_idx=hook['peak_idx'] if hook['hook_effect'] else None)
pd.DataFrame([fit]).to_csv('fit.csv', index=False)
```

The helper is `scripts/cooperativity_dc50.py`. Inputs are measured concentrations and percent degradation; concentrations must be in nM for the helper's fit bounds and caveat text. The source hook flag uses a 15-percentage-point downturn and the fit caveats insufficient plateau. Alpha uses matched measured binary and ternary Kd values. See `references/measurements.md`.

Done when assay conditions, fitted values, hook status and caveats are reported.
