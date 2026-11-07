# Measure gene and site concordance
~~~bash
iqtree3 -t concat.treefile --gcf loci.treefile --prefix concord_g
~~~
Use one compatible gene tree per locus and a reference tree. Separately compute likelihood site concordance:
~~~bash
iqtree3 -te concat.treefile -s concat.fasta --scfl 100 -T 4 --prefix concord_s
~~~
The source uses separate calls because the flags cannot be combined; --scfl requires a compatible version. Low gene concordance with strong bootstrap indicates discordance, not by itself proof of ILS or introgression.
See [support and concordance](../references/support.md).
Done when decisive loci/sites, discordance and support are interpreted together.
