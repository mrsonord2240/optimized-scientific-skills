# Audit a palette

```bash
Rscript scripts/audit_palette.R colors.csv palette-audit.csv
```

Input CSV must contain a `color` column of valid hex colors; supply the actual ordered mapping. The output reports minimum pairwise CIELAB distances under deutan, protan and tritan simulations plus a luminance monotonicity verdict. It overwrites the output file. Inspect the actual figure and simulations too: a minimum distance is not an accessibility guarantee and grayscale and CVD checks address different issues. Categorical colors may need shapes, direct labels or facets; preserve names across subsets. Document palette and type, normalization and reference, clipping, missing count and encoding, distances, luminance, redundant encoding and visual inspection limitations.

Done when the actual mapping and figure have been inspected and the audit record includes these details.
