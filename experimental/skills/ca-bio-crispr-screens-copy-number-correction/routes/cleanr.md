# Apply CRISPRcleanR
~~~bash
Rscript scripts/run_crispr_cleanr.R
~~~
Run from a directory containing counts.txt with sgRNA, gene and sample columns. The source example uses built-in KY_Library_v1.0; adapt the script to actual guide coordinates and count format before use. It uses min_reads=30, saves an RData result, and optionally reads cn_profile.txt with gene and copy_number. Its commented downstream filename is not generated: export and verify a TSV explicitly before downstream use. Example sample variables are not passed to normalization. Check package version and inspect correction and residual focal bias. See [method details](../references/methods.md).
Done when corrected results are exported, matched to inputs, and validated against uncorrected results and known controls.
