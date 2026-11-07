# Conformal regression baseline

```python
import numpy as np
from mapie.regression import MapieRegressor
from sklearn.ensemble import RandomForestRegressor
X_train, y_train, X_test = np.load('X_train.npy'), np.load('y_train.npy'), np.load('X_test.npy')
model = MapieRegressor(estimator=RandomForestRegressor(n_estimators=500, random_state=42, n_jobs=-1), method='plus', cv=5)
model.fit(X_train, y_train)
pred, intervals = model.predict(X_test, alpha=0.1)
np.savez('conformal.npz', predictions=pred, intervals=intervals)
```

This source example targets MAPIE >=0.8,<1.0 and an sklearn baseline, not a chemprop wrapper. Nominal 90% coverage depends on assumptions; assess empirical coverage on a deployment-relevant holdout. See `references/uncertainty.md`.

Done when predictions, intervals and empirical coverage are reported.
