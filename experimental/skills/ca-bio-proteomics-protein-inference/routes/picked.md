# Estimate picked group FDR
~~~bash
python scripts/picked_group_fdr.py groups.tsv --decoy-prefix DECOY_ --score-col probability --out passing.tsv
~~~
Input needs semicolon-separated accessions and higher-is-better scores. The source helper pairs exact accession sets; mismatched target/decoy group memberships remain unpaired, and fewer than ten decoy groups trigger a warning. Its threshold is fixed at 0.01. Inspect decoy prefix, pairs, group resolution and output; see [picked FDR](../references/picked-group-fdr.md).
Done when target/decoy coverage, pair behavior and group-level results are reviewed.
