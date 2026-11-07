# Compare candidate topologies
~~~bash
iqtree3 -s alignment.fasta -m GTR+G -z trees.nex -n 0 -zb 10000 -au --prefix autest
~~~
Supply the reviewed model and a priori candidate trees; the source uses 10000 RELL replicates. An AU p-value below 0.05 rejects a candidate under the test; failure to reject is not acceptance. Confirm current version and output.
Done when candidate set, model and test uncertainty are reported.
