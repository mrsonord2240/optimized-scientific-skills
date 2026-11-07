# Single-parameter FDR-aware sample size
~~~sh
Rscript scripts/single_size.R scenario.csv result.rds
~~~
Supply a one-row CSV with nGenes, pi0, m, mu, disp, fc, fdr, power, maxN and seed. Here m is the source's internal pseudo-sample size, not replicates per group. The script uses ssizeRNA_single, preserves its full result and stops if no n within maxN is returned. These are model-dependent estimates; compare plausible effects and dispersion and inspect realized FDR under the intended analysis. Verify installed ssizeRNA APIs.
**Done when:** the seed, inputs, returned size and uncertainty are documented, including unreachable scenarios.
