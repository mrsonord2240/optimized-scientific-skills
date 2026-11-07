# Pinned-release ORA

```bash
Rscript scripts/wikipathways_pinned_enrich.R sig_entrez.txt universe_entrez.txt 'Homo sapiens' wp_pinned.csv
```

Run from the skill directory with one Entrez ID per line in each input; the universe must represent all tested genes. The organism is the exact WikiPathways scientific name and the output argument is optional. The script contacts the WikiPathways archive, tries recent dated releases, prints the date selected and writes results if a CSV path is supplied. Record that date; it selects a recent available release, not a user-specified historical release. Archive availability and package interfaces need verification at use time. If no candidate release is available, stop and obtain an appropriate archived GMT separately; the script does not implement that fallback. Check mapping and empty results before interpreting them.

**Done when:** the output and reported release date are retained with the tested universe and mapping provenance.
