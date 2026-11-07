# Exploratory penalized Cox fit
~~~sh
python scripts/cox_fit.py training.csv time event coefficients.csv --alpha 0.1 --l1-ratio 0.9
~~~
CSV rows are independent subjects; time is nonnegative follow-up, event is 0/1, and remaining columns are complete numeric predictors prepared using training data only. This fits a fixed-penalty research model and writes coefficients; it does not tune parameters, estimate performance or produce a validated predictor. Preprocessing, feature selection and penalty tuning must occur within training folds for any held-out evaluation. Verify proportional hazards, censoring and installed API.
**Done when:** fit assumptions, data provenance, penalty and independent validation plan are documented.
