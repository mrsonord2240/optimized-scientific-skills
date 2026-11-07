# Other study designs

```bash
Rscript examples/aldex2_analysis.R
```

This source example uses synthetic counts and writes a fixed `da_consensus.csv` in the working directory; use a fresh directory. The real-data wrapper is limited to two independent groups. ANCOM-BC2, LinDA mixed models, MaAsLin2/3, ZicoSeq and LEfSe need input- and version-specific adaptation. For ANCOM-BC2 the source sets BH explicitly rather than the reported Holm default and checks pseudo-count sensitivity; model coefficient names vary. Match random effects to actual repeated subjects and cross-check implementations. LEfSe is exploratory and not an FDR-controlled consensus substitute. Declare and test prevalence filters. See [assumptions](../references/assumptions.md).

Done when the model matches the experimental unit and all methods, coefficients, correction families and disagreements are documented.
