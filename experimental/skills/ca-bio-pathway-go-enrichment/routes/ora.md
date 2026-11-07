# Gene list plus universe

```bash
Rscript scripts/ora_go.R genes.txt universe.txt results.csv keytype=SYMBOL ont=BP
```

`genes.txt` and `universe.txt`: one ID per line. `universe.txt` holds every gene that could have entered the list (expressed or tested genes, detected proteins, panel genes). The script maps both with the same call, drops duplicates, runs `enrichGO` (BH, `minGSSize=10`, `maxGSSize=500`, `pvalueCutoff=0.05`, `qvalueCutoff=0.2`) and runs `simplify` (cutoff 0.7) per ontology.

- `keytype=` is one of `SYMBOL`, `ENTREZID`, `ENSEMBL` or any key of the OrgDb. Other species: `orgdb=org.Mm.eg.db`, `org.Dr.eg.db`, `org.Sc.sgd.db` (with `keytype=ORF`).
- `ont=ALL` runs BP, MF and CC separately and simplifies each; the output has an `ONTOLOGY` column.
- `simplify=false` keeps every term.
- Read `results.csv.log`: report the foreground and universe sizes and the conversion rate. Over 15% of genes lost makes the result unreliable.
- A list under about 15-20 genes has little power: report counts and fold enrichment and consider GSEA (`routes/other-skills.md`).

Done when `results.csv` exists and the answer reports `p.adjust`, `FoldEnrichment`, `Count`, and the universe size. If it is empty, see `routes/errors.md`.
