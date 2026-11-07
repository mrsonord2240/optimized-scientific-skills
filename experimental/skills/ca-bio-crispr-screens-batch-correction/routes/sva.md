# Sva
~~~bash
Rscript scripts/sva.R counts.tsv metadata.tsv design.tsv
~~~
Use feature-indexed counts and sample-indexed metadata containing condition. The full model protects condition and the null model has an intercept. Review factors and design rank; surrogate variables cannot resolve full confounding. See [SVA notes](../references/sva.md).
Done when the design matrix is saved and its suitability checked.
