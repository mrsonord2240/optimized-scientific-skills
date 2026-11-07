## Conformal Prediction for Calibrated Uncertainty

Use conformal prediction when calibrated marginal coverage under the stated exchangeability assumptions matters. Ensemble variance is simpler, but it is not a substitute for a conformal guarantee.

```python
# MAPIE expects a scikit-learn-compatible estimator (.fit / .predict / .predict_proba).
# chemprop 2.x is NOT scikit-learn-compatible out of the box -- either wrap chemprop
# in a thin sklearn estimator class or use MAPIE only with the sklearn baseline.
from mapie.regression import MapieRegressor
from sklearn.ensemble import RandomForestRegressor

base = RandomForestRegressor(n_estimators=500, random_state=42, n_jobs=-1)
mapie = MapieRegressor(estimator=base, method='plus', cv=5)
mapie.fit(X_train, y_train)
y_pred, y_intervals = mapie.predict(X_test, alpha=0.1)  # alpha=0.1 -> 90% coverage
```

Alpha 0.05 targets 95% marginal coverage and alpha 0.10 targets 90%, subject to the conformal method's assumptions. MAPIE supports the sklearn baseline directly; integrating chemprop requires a separately implemented and tested compatible wrapper.

`method='plus', cv=5` fits the base estimator six times, so this is expensive for large fingerprint matrices even with `n_jobs=-1`. Use a cheaper split-conformal design when repeated base-model fitting is not affordable, and always report held-out empirical coverage for the selected split.

## SHAP / Atomic Attribution

For mechanistic interpretation:

For a scikit-learn-style model (e.g., Random Forest baseline on ECFP4), SHAP integrates directly:

```python
import shap
from sklearn.ensemble import RandomForestClassifier

# X_train / X_test are Morgan fingerprint arrays (n_samples, n_bits)
model = RandomForestClassifier(n_estimators=500, random_state=42).fit(X_train, y_train)
explainer = shap.TreeExplainer(model)
shap_values = explainer.shap_values(X_test)

# Per-bit contribution; for atomic interpretation, map bits back to
# generating atoms via AllChem.GetMorganFingerprintAsBitVect(mol, ..., bitInfo=bi)
# and aggregate SHAP across all bits triggered by each atom.
```

For chemprop D-MPNN, SHAP requires a custom wrapper (chemprop is not sklearn-compatible). A PyTorch attribution method must be adapted to the model's graph inputs and validated; the chemprop 2.x CLI does not provide the `--uncertainty-method classification` atom-attribution interface. Use directly supported fingerprint SHAP for the classical baseline unless a tested graph-attribution implementation is available.

## Bayesian Optimization for Active Learning

```python
import numpy as np
from scipy.stats import norm
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF

gp = GaussianProcessRegressor(kernel=RBF(length_scale=1.0), random_state=42)
gp.fit(X_train, y_train)
mu, sigma = gp.predict(X_pool, return_std=True)

# Expected Improvement
def expected_improvement(mu, sigma, y_best, xi=0.01):
    improvement = mu - y_best - xi
    ei = np.zeros_like(mu, dtype=float)
    nonzero = sigma > 0
    z = improvement[nonzero] / sigma[nonzero]
    ei[nonzero] = (
        improvement[nonzero] * norm.cdf(z)
        + sigma[nonzero] * norm.pdf(z)
    )
    return ei

ei = expected_improvement(mu, sigma, y_train.max())
next_to_test = X_pool[ei.argmax()]
```

For chemprop + active learning, replace GP with chemprop ensemble + ensemble variance.

## Calibration (Platt / Isotonic)

Deep learning probabilities are not guaranteed to be calibrated. Use Platt (logistic) or isotonic calibration for binary probabilities, choosing the method with a held-out calibration set. `--metric roc` evaluates ranking and does not automatically calibrate chemprop probabilities; export validation probabilities and fit the calibrator externally:

```python
from sklearn.isotonic import IsotonicRegression
iso = IsotonicRegression(out_of_bounds='clip').fit(val_chemprop_probs, val_true)
test_calibrated = iso.predict(test_chemprop_probs)
```
