# Inspect scalar convergence
~~~bash
python scripts/bayesian_convergence.py run1.p run2.p 0.25
~~~
Use two independent MrBayes parameter files with matching headers. The source helper uses an approximate ESS and PSRF and 25% illustrative burn-in; compare with established diagnostics and inspect traces. The helper refuses PSRF when within-chain variance is zero; inspect constant parameters explicitly. Passing scalar checks does not establish topology convergence: assess ASDSF, max split discrepancy or RWTY separately. See [diagnostic notes](../references/methods.md).
Done when all parameters, topology, burn-in and prior sensitivity have been reviewed.
