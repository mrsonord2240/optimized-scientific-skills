# KEGG over-representation

```bash
Rscript scripts/kegg_ora.R hits.txt universe.txt ORGANISM KEYTYPE result.csv
```

Each input text file has one ID per line in the same KEGG-supported namespace. Replace organism and key type after verifying the installed API and the organism/strain. The wrapper requires a nonempty measured universe and hits within it; it uses the source settings BH, p cutoff 0.05, q cutoff 0.2 and set sizes 10–500. `enrichKEGG` queries the live KEGG service; record access date, mapping loss and database version. Output may be overwritten. For a reproducible snapshot, the source GSON recipe requires adapting the retained `examples/kegg_enrichment.R` and retaining the snapshot; no snapshot is created by the wrapper.

Done when ID mapping, eligible universe, tested pathways, reference date and adjusted results are reviewed.
