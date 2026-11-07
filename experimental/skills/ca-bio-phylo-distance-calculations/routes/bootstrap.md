# Bootstrap rooted identity NJ trees
~~~bash
python scripts/bootstrap_consensus.py alignment.fasta OUTGROUP consensus.nwk --replicates 50
~~~
Replace OUTGROUP with an exact taxon name. Source uses identity distances, seed 42 and majority cutoff 0.5; this wrapper does not strip gaps, so prefilter alignment consistently before use. Every replicate is rooted on the same outgroup to avoid splitting support between complementary rooted clades.
Done when outgroup, gap policy, model and support limitations are reviewed.
