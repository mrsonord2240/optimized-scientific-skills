# Fit multi-condition MLE
~~~bash
mageck mle --count-table screen.count.txt --design-matrix design.txt --output-prefix timecourse_mle --norm-method median --permutation-round 10
~~~
Supply a sample-matched design with baseline and condition or covariate columns. Check full rank, replication and confounding; adding batch cannot solve full confounding. For optional guide efficiencies, see [model notes](../references/models.md). NaN beta is undetermined, not zero. The source reports default two permutation rounds and recommends 10 for borderline calls; report the setting, assess stability and compare wald-fdr without claiming either validates results.
Done when matrix, effects, missing estimates and near-threshold stability are reviewed.
