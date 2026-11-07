# Reactome ORA

```bash
Rscript scripts/reactome_ora.R hits.txt universe.txt human results.csv
```

Inputs contain one Entrez ID per line, in the same namespace, with hits a subset of the measured eligible universe. Replace organism with one supported by the installed ReactomePA and its annotation package. The wrapper uses BH, p cutoff 0.05, q cutoff 0.2 and set sizes 10–500, and reports a NULL result rather than treating it as a successful enrichment. It uses installed local annotation; record package versions and mapping losses. Output may be overwritten. Parent and child pathways share genes and should not be counted as independent findings.

Done when mappings, universe, annotation version, tested results and hierarchy overlap are reviewed.
