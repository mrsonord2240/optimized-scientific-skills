# Cox and MR displays

```bash
Rscript examples/forest_phd.R
```

This source demonstration uses synthetic study and Cox data and package-supplied MR data; it writes fixed files under `plots/` and is not an analysis of user data. Adapt its inputs and model to the actual study before use. Adjusted covariate Cox plots are not subgroup-treatment forests. For prespecified subgroups use stratum estimates and a formal treatment-by-subgroup interaction test; the source stops for fewer than 20 observations or five events in a stratum. For MR compare methods, print the MR-Egger intercept, and do not treat agreement as proof of causality. A generic input-file interface for these model-specific alternatives is not supplied. Source also describes forestplot table layouts, netmeta and chronological cumulative meta-analysis; they require input-specific adaptation.

Done when the intended model, group sizes, interaction or method diagnostics and interpretation limitations are documented.
