# Evaluate Chronos and CERES
**Source-incomplete:** The source supplies a Chronos fragment with undefined input dataframes, not a complete file-based executable pipeline; CERES is described without a runnable recipe. Supply and validate required data and adapt [panel methods](../references/panel.md).
Chronos training uses sample-by-guide readcounts, sequence map, guide-to-gene map and nonempty negative controls; the example uses 301 epochs and seed 20260924. Post-training alternate_CN requires at least three cell lines and a matched line-by-gene CN matrix. Validate each line and amplicon; a panel average can hide residual bias.
Done when input schemas, version-specific API and corrected gene effects are independently validated, or the method is explicitly deferred.
