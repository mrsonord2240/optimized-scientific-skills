# Assess concordance
~~~bash
iqtree3 -t species.tre --gcf gene_trees.nwk --prefix cf_g
~~~
For likelihood site concordance, use a separate call:
~~~bash
iqtree3 -te species.tre -s concat.fasta --scfl 100 --prefix cf_s
~~~
Use compatible taxa and model/version. localPP is not bootstrap; gCF and sCF represent different data summaries. Minority-quartet asymmetry can motivate further tests but is not alone proof of introgression or its absence. Inspect gene-tree uncertainty and other model violations.
Done when discordance, decisive data and uncertainty are reported.
