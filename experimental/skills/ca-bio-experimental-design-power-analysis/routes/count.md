# Negative-binomial approximation
~~~sh
Rscript scripts/count_power.R scenarios.csv power.csv
~~~
Each CSV row needs depth, cv, effect, n, power and alpha columns; leave exactly one of effect, n or power empty for RNASeqPower to solve. Depth is per-feature coverage, not total library reads. Values must be finite and in the stated domains; n returned is per group and must be rounded up for planning. This is a single-feature approximation, not transcriptome-wide FDR power. Check package version and report plausible sensitivity scenarios.
**Done when:** results, inputs, units and assumptions are documented and checked against an appropriate study-specific simulation.
