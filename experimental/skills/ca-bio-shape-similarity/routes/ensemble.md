# Conformer-ensemble shape search

```bash
python scripts/shape_search_ensemble.py --query "$(head -n 1 query.smi | awk '{print $1}')" --library library.smi --n-conf 20 > shape_results.txt
```

For an alternative pairwise multi-conformer helper see `scripts/shape_search.py`. The script reports dropped structures and nonconverged conformers in its output. Add --scaffold-hop with calibrated --shape-threshold and --ecfp-threshold for shape-high/2D-low candidates. It embeds one query conformer and an ensemble for each target; assess coverage. See `references/open3dalign.md` and `references/interpretation.md`.

Done when scores, warnings and search settings are saved and reviewed.
