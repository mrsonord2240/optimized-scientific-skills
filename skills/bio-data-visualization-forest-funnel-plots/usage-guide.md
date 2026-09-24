# Forest and Funnel Plots - Usage Guide

Ask for the analysis unit and the output question, for example:

- "Pool these 15 study-level log-ORs with REML and show a log-scale forest, Q, tau-squared, I-squared, and a prediction interval."
- "Make a contour-enhanced funnel; run Egger only if there are at least 10 studies, and treat trim-and-fill as sensitivity-only."
- "Show treatment HRs by prespecified subgroup, test the treatment-by-subgroup interaction, and stop on sparse strata."
- "Compare IVW, weighted median, mode, and MR-Egger in an MR forest without SNP rows; print the MR-Egger intercept."

Provide effect estimates and standard errors (or enough raw data to derive them), the effect scale, study labels, and any prespecified subgroup definition. A forest of adjusted Cox covariates is not a subgroup-treatment forest; request the interaction analysis when effect modification is the question.

The runnable reference is [examples/forest_phd.R](examples/forest_phd.R). It uses labelled synthetic data for the meta-analysis and Cox demonstration, then the package's bundled lipid/CHD data for MR. Replace those inputs with your own data before drawing substantive conclusions.
