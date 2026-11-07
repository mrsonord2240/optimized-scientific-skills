# Summarize observations by experimental unit
~~~sh
Rscript scripts/summarize_units.R observations.csv unit_means.csv
~~~
Input CSV requires unit, condition and numeric measurement columns. The helper checks that each unit has only one condition and reports the mean and observation count by unit. Verify that unit corresponds to independent treatment assignment; aggregation is not a substitute for a justified hierarchical or repeated-measures model in every design.
**Done when:** experimental-unit counts by condition and the aggregation choice are documented.
