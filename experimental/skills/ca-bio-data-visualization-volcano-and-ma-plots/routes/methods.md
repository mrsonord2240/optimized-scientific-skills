# Shrinkage and alternatives

```bash
Rscript scripts/volcano_csv.R results.csv volcano.pdf
```

This command plots an existing table and does not perform shrinkage. The source `scripts/volcano_phd.R` is a session-dependent recipe requiring a fitted `dds` and the named coefficient; adapt those inputs and identifiers before use. It demonstrates DESeq2 apeglm coefficient shrinkage, ashr contrast shrinkage, MA and EnhancedVolcano. An `svalue=TRUE` result has s-values instead of p-value and adjusted-p columns; retain the corresponding original results if both are needed. EnhancedVolcano's y quantity and threshold must match; verify requested labels against identifiers and missing values. These alternatives are not standalone input-file workflows here. See [interpretation](../references/interpretation.md).

Done when the chosen estimate, comparison, significance quantity and any shrinkage settings are documented and reviewed.
