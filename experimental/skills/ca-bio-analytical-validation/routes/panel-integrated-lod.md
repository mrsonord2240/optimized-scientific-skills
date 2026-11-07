# Multi-locus sampling lower bound

```bash
python scripts/panel_integrated_lod.py --input-ng 30 --n-loci 16 32 48 --min-loci 2 --out panel_bound_results.tsv
```

For each panel size, the lowest VAF on a log grid (1e-6 to 1e-2) where the chance that at least `--min-loci` loci have a template present reaches `--target` (default 0.95). Output columns carry `vaf_sampling_lower_bound` and the interpretation string. The script prints the four assumptions: equal VAF at every locus, independent sampling, perfect recovery and detection of present templates, no depth failures, background error, false positives or calling effects. It refuses `--min-loci` above `--n-loci`, and a grid that never reaches the target.

- This is an optimistic best case. Never call it LoD95. Panel LoD95 is established empirically on replicated full-workflow dilutions with the k-of-N rule fixed before fitting (`routes/lod95-probit.md`).
- Multi-locus integration can improve sensitivity, but real loci differ in recovery and background, and detections are not independent.
- Single hotspot versus bespoke panel at low VAF: use this only to explore the sampling bound, then validate the chosen panel and rule.
- Input mass is in `--input-ng`; `--ge-per-ng` overrides 303.03 and must be stated.

Done when `panel_bound_results.tsv` exists and the answer says "theoretical sampling lower bound" with the assumptions.
