# Run a configured MrBayes analysis
~~~bash
mb analysis.nex
~~~
Supply an aligned NEXUS containing a reviewed MrBayes block. The source example uses GTR+I+G, a compound-Dirichlet branch prior, two runs, four chains, 10 million generations and 25% illustrative burn-in; stoprule=no avoids halting solely on topology. Settings must fit the data and available tools; see [configuration](../references/methods.md). For an existing BEAST2 XML the source command is beast -seed 12345 -overwrite my_analysis.xml; choose a fresh output to avoid overwriting results. RevBayes scripts require a bespoke model.
Done when settings, priors, outputs and independent convergence diagnostics are documented.
