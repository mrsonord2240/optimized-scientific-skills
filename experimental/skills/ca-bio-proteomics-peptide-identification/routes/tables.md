# Estimate PSM FDR from tables
~~~bash
python scripts/table_fdr.py search_results.tsv --scan scan --score score --protein protein --out psms_1pct.tsv
~~~
Use this for concatenated target-decoy competition, a higher-is-better score and unique spectrum identifiers across all runs. It keeps the highest hit per scan and uses (decoys+1)/targets; verify recognized prefixes and sufficient decoys. For separate target and decoy searches use:
~~~bash
python scripts/separate_search_fdr.py target.tsv decoy.tsv
~~~
The source separate-search example reports counts for estimated pi0 and pi0=1; it does not export a filtered file. Use its functions in a custom workflow if a table is required. See [FDR methods](../references/fdr_from_tables.md). q-values characterize list-level error; PEP is a distinct per-match quantity.
Done when search mode, decoy ratio, score direction and confidence level are documented.
