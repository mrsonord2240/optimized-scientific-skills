# Nested cross-validation
~~~sh
python scripts/nested_cv.py samples.csv outcome cv_results.json
~~~
CSV rows must be independent samples; the named binary outcome column is excluded from numeric features. The source pipeline standardizes and selects features inside each fold, with nested tuning of k. This wrapper uses stratified folds and is not appropriate for repeated donors, families, sites or time dependence; use group-aware or temporal splits when required. Define preprocessing, feature selection and the model entirely within training folds, and evaluate an independent cohort.
**Done when:** fold metrics, split design, tuning and uncertainty are reported without treating them as external validation.
