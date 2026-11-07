# Matrix and per-sample alternatives

```bash
Rscript examples/gsea_go.R
```

This source example uses a synthetic ranking and local annotation and writes a temporary result; it is not a user-data matrix test. CAMERA, fry, GSVA and ssGSEA require a matching expression matrix, design and gene sets and need input-specific adaptation. CAMERA and fry test different hypotheses; GSVA/ssGSEA produce per-sample scores rather than automatic pathway contrast p-values. Source recommends assessing inter-gene correlation and setting GSVA's kernel for the input scale. The retained MSigDB example also uses synthetic data and version-specific package fields. Verify APIs and clean process exit.

Done when model design, gene-set mapping, null, score interpretation and assumptions are documented for the actual data.
