# Simulation and donor-level power
For the source synthetic bulk example, run:
~~~sh
Rscript scripts/simulation_example.R
~~~
For the synthetic donor demonstration, run:
~~~sh
Rscript scripts/pseudobulk_example.R
~~~
Both use fixed illustrative parameters and are not study-specific estimates. The source PROPER workflow assumes illustrative empirical priors and fixed simulation settings; a study-specific pilot mean-dispersion trend, plausible effect distribution and intended test must be supplied and validated before using it as a design estimate. Its pseudobulk example is synthetic and uses donor-level replication; adapt it to the study's donor and cell structure, not cells as independent patients. The source proposes checking simulated realized FDR against a prespecified target and tolerance before selecting power results; account for Monte Carlo uncertainty, particularly with few simulations. See [simulation guidance](../references/methods.md). The source does not provide a general file-based pilot-model interface for PROPER or powsimR.
**Done when:** assumptions, simulation settings, uncertainty and realized error rates are reviewed before selecting a design.
