# Background

Not needed to run a route.

In p much larger than n, a flexible learner's variance dominates, full covariance is singular so QDA and full LDA are undefined, and simple diagonal or linear methods match or beat elaborate ones. Random forest is not the obvious choice for expression data; SVM and regularized logistic frequently win on microarray-style problems, and gradient boosting beats deep nets on tabular omics. Regularization is load-bearing.

In diagnostic use the probability is the product, so calibration breaks silently: a model can rank well (AUC 0.9) and still output dishonest risks. The commonest cause of a beautiful AUC is a batch artifact: when batch correlates with outcome the model learns the cleaner technical signal.

Measured figures behind the routes (synthetic data with known generators unless Golub is named): batch confounding with zero signal gave random-split AUC 0.82, 0.88, 0.86 for 2, 3, 6 batches against leave-one-batch-out per-batch means 0.54, 0.70, 0.45, batch predictability 1.00. A 1,500-sample, 8% prevalence logistic: unweighted predicted 0.99x observed (slope 0.93, AUC 0.77); balanced 4.2x (slope 0.72); balanced then recalibrated on 1,000 held-out samples 1.11x. Random forest AUC plain 0.776 against balanced 0.812 (higher in 6 of 6 seeds); logistic 0.826 against 0.820. XGBoost at 300 samples stopped at round 0 (AUC 0.60); at 600 samples reached 0.83.

Failure modes: beautiful AUC that is batch; forest or boosting probabilities trusted as risks; class weights or SMOTE for a risk model, or SMOTE before the split; early-stopping or calibration data reported as the result; out-of-bag error after leaky selection.
