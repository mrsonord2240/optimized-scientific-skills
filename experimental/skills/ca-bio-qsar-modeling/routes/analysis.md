# Fingerprint attribution

```python
import numpy as np, shap
from sklearn.ensemble import RandomForestClassifier
X_train, y_train, X_test = np.load('X_train.npy'), np.load('y_train.npy'), np.load('X_test.npy')
model = RandomForestClassifier(n_estimators=500, random_state=42).fit(X_train, y_train)
values = shap.TreeExplainer(model).shap_values(X_test)
np.save('shap_values.npy', np.asarray(values))
```

Inputs must use the same fingerprint representation and split. Bit-to-atom mapping requires separate interpretation; chemprop graph attribution requires a validated custom implementation. The source Gaussian-process expected-improvement example is in `references/uncertainty.md`.

Done when attribution and its representation limits are recorded.
