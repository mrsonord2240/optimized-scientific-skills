# Infer Bayesian group posteriors
~~~bash
python scripts/epifany_inference.py peptides_with_pep.idXML --out groups.tsv
~~~
Input PSMs need posterior error probabilities and target/decoy evidence. Output is group posterior probabilities, not a filtered FDR result; use [picked groups](picked.md) only after checking group resolution. Source notes --greedy-group-resolution did not remove subsumable groups in its example; inspect them. See [EPIFANY notes](../references/epifany-bayesian-inference.md).
Done when posterior input, groups and resolution limitations are reviewed.
