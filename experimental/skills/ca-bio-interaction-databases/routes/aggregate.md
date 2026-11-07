# Aggregate interactions

```bash
python scripts/query.py --resource aggregate --genes genes.txt --out aggregate.tsv
```

The source aggregate is an undirected union among query genes of STRING at 700, OmniPath post-translational interactions and optional BioGRID low-throughput physical evidence when `BIOGRID_API_KEY` is set. It excludes self-interactions, collapses direction and uses OmniPath's academic default; do not use it as a signed or commercially cleared graph. Overlap across resources is not calibrated independent confidence. For signed direction use the separate curated tables; see `references/resources.md`.

`scripts/interaction_clients.py` contains the reusable source clients and graph helpers; the aggregate uses its default STRING caller identity.

Done when the edge table and per-edge sources are saved and scope and direction loss are documented.
