# Compare models by stepping-stone
~~~bash
mb model_comparison.nex
~~~
Supply a NEXUS file containing appropriate model, proper priors and the source stepping-stone command ss ngen=1000000 nsteps=50 samplefreq=100 diagnfreq=1000. Run each comparable model separately and record marginal log likelihoods; 2 ln BF is twice their difference. Check sampling and precision; do not select models from the harmonic-mean estimator.
See [model comparison and CAT](../references/model-comparison.md).
Done when prior choices, sampling and model comparison uncertainty are documented.
