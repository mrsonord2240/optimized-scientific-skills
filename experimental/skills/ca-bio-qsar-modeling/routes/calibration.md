# Probability calibration

```python
import pandas as pd
from sklearn.isotonic import IsotonicRegression
validation = pd.read_csv('validation.csv')
test = pd.read_csv('test.csv')
model = IsotonicRegression(out_of_bounds='clip').fit(validation['probability'], validation['label'])
test['calibrated'] = model.predict(test['probability'])
test.to_csv('calibrated.csv', index=False)
```

Use a separate calibration set and report calibrated and uncalibrated performance. Source also discusses Platt scaling; see `references/uncertainty.md` and `references/errors.md`.

Done when calibration and held-out metrics are reported.
