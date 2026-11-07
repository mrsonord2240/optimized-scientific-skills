# FUSION association and conditional analysis
~~~sh
/path/to/skill/scripts/fusion_twas.sh gwas.sumstats weights.pos weights_dir/ ld/EUR. out/twas "22"
~~~
The script is scripts/fusion_twas.sh. Run from the FUSION clone with an existing output directory. Inputs need SNP A1 A2 Z, matching weight files and ancestry-matched chromosome LD reference. The optional sixth argument is a space-separated chromosome list; it defaults to 1-22. The script runs association per chromosome and post-processing with a 100000 bp window, overridable via LOCUS_WIN. Check installed CLI and the source-reported top1 single-SNP post-processing issue before running. [Pipeline details](../references/fusion-pipeline.md); [failure modes](../references/failure-modes.md).
**Done when:** chromosome and combined association tables and conditional results are produced and checked; conditional Z values are not FOCUS PIPs.
