# STRING network

```bash
python scripts/query.py --resource string --genes genes.txt --out edges.tsv --caller research-project --threshold 700 --network-type functional
```

Supply one gene symbol per line and the correct taxon (`--taxon`, default 9606). Functional combined scores include nonphysical evidence; select `--network-type physical` for STRING's physical-complex network. The experiment channel alone is not a physical-only filter. For identifiers use `--resource resolve`. Record STRING version, thresholds and evidence channels; see `references/resources.md`.

Done when the edge or ID table is saved and scope, species, version and evidence are reported.
