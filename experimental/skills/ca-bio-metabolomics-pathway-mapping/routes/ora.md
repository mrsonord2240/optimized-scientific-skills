# Identified-compound ORA

```bash
Rscript scripts/local_ora.R hits.txt assay_universe.txt ora.csv 2
```

Inputs are one KEGG compound ID per line for hits and the assay-coverage universe; check that hits belong to the universe. The script fetches a public KEGG reference table over the network, computes hypergeometric values locally and applies BH across all nonempty pathways before filtering the output to at least two hits. Its source converts generic `map` pathway identifiers to `hsa` labels; this does not establish organism-specific membership and requires independent verification. Mapping with `Rscript scripts/map_compounds.R compounds.txt hsa name kegg_ids.txt` may download public reference libraries; inspect mapping coverage. The alternative `scripts/ora_api.R` transmits user compound data to a remote service and the source observed filtered requests can fail. Obtain appropriate authorization before using it. Files can be overwritten.

Done when mapping coverage, actual pathway reference, universe, full testing family and annotation uncertainty are reported.
