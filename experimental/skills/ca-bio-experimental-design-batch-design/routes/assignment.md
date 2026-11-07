# Constrained batch assignment
~~~sh
Rscript scripts/assign_batches.R samples.csv 3 8 condition,sex layout.csv 10000 1
~~~
CSV needs id and the named feature columns. Positional arguments specify batches, positions per batch, comma-separated features, output, optional iterations (default 10000) and seed (default 1). The source designit API is version-sensitive; verify the installed version. Its coverage check intentionally stops if any level is absent from any batch, a conservative requirement stronger than estimability. An incomplete-block design may still be identifiable; assess the full model matrix rank and design constraints rather than calling every empty cell perfect confounding. Inspect balance and randomization after optimization.
**Done when:** the layout satisfies capacity and study constraints, balance tables and model estimability have been reviewed, and seed and metadata are recorded.
