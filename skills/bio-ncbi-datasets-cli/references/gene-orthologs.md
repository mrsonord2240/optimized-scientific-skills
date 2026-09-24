# Gene records across species and ortholog sets

Read when the request spans more than one species or asks for orthologs.

### Gene metadata across species

`--taxon` on `summary gene symbol` is **single-species only** (it picks which species' gene record
to resolve the symbol against; default `human`) -- it does not accept a clade like `Mammalia` and
errors outright if you try (`gene requires an at-or-below-species-level taxon`). The only mechanism
this subcommand has for a genuinely multi-species pull is `--ortholog <taxon|all>`, which accepts
any taxonomic rank (not just `all`) and returns NCBI's ortholog set for that clade -- one
representative gene per species, limited to vertebrates and insects:

```bash
datasets summary gene symbol BRCA1 \
    --ortholog Mammalia \
    --as-json-lines \
  | dataformat tsv gene --fields gene-id,symbol,tax-name,description,chromosomes \
  > brca1_mammals.tsv

head brca1_mammals.tsv
```

Verified live (18.37.0): this returns 272 real rows across Mammalia (human, mouse, rat, dog, cow,
macaque, chimp, opossum, pig, ...). Outside vertebrates/insects, or for a single specific species,
loop `--taxon <species>` per species instead.

### Find orthologs for a gene

```bash
datasets summary gene symbol BRCA1 --taxon human --ortholog all --as-json-lines \
  | dataformat tsv gene --fields gene-id,symbol,tax-name,description \
  > brca1_orthologs.tsv
```

`--ortholog` takes a required value (`all`, or one or more taxa) -- a bare `--ortholog` flag is
consumed as swallowing the next flag's value and fails with a misleading "taxonomy name not exact"
error. It returns NCBI's ortholog set (a single representative per species; tree-aware orthology
with multiple co-orthologs is in `ortholog-inference` / Compara / OMA).
