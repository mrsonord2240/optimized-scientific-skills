<!-- Moved from SKILL.md (2026-09-21); code now in scripts/. -->

## Network-Diffusion Enrichment (FELLA)

**Goal:** Return the intermediate enzymes, reactions, and modules that mechanistically link the affected metabolites, not just a ranked pathway list.

**Approach:** Build the KEGG knowledge graph once, then per-analysis map KEGG IDs and run heat diffusion; inspect excluded (unmapped) compounds explicitly.

```bash
# cpd_ids.txt: KEGG compound IDs only, one per line; builds the fella_hsa database on first use, reuses it after
Rscript scripts/fella_diffusion.R cpd_ids.txt fella_hsa fella_results.csv
```

The graph build hits the live KEGG API (slow), so the database is cached in the directory. Unmapped compounds are printed via `getExcluded()`: report them. `diffusion` is the recommended default; `runHypergeom` = plain ORA over the graph, `runPagerank` (lowercase r) = directed random walks; the method string is lowercase. `approx = 'normality'` (used) is analytic/deterministic, no seed needed; with `approx = 'simulation'` call `set.seed()` first, since it resamples `niter` times.
