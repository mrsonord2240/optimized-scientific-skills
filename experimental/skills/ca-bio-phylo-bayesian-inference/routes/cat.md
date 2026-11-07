# Diagnose PhyloBayes chains
~~~bash
bpcomp -x 1000 10 chain1 chain2
~~~
Use independent completed chain files and a justified burn-in/thinning interval. For scalar traces:
~~~bash
tracecomp -x 1000 chain1 chain2
~~~
The source reports maxdiff <0.1 and effsize >300 as good conventions; inspect all diagnostics and chains, not just thresholds. CAT-GTR run settings and caveats are in [methods](../references/model-comparison.md).
Done when scalar and topology agreement and retained sample choices are reviewed.
