# Curated and signed interactions

```bash
python scripts/query.py --resource signor --genes genes.txt --out signor.tsv
```

SIGNOR records are directed and signed with mechanism where available; preserve those fields and report repeated site or paper records. For BioGRID low-throughput physical records select `--resource biogrid` and set `BIOGRID_API_KEY` in the environment; do not print the key. For OmniPath use `--resource omnipath --license academic` or `commercial`; its source-level commercial screen is not a substitute for reviewing remaining terms and derived fields. See `references/resources.md` and `references/terms.md`.

Done when the requested interaction table is saved with source, evidence, direction and applicable usage terms.
