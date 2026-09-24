<!-- Moved from SKILL.md "ORA on an Identified Compound List" (2026-09-21); code now in scripts/. Run scripts/map_compounds.R (SKILL.md) first: it produces `kegg_ids.txt`. -->

### Local-Only ORA (default: no remote call on user data, verified background correction)

Computes the hypergeometric ORA locally from KEGGREST's public pathway-to-compound table (generic
reference data, not user data). Checked on KEGGREST 1.46.0: ~3 s to fetch and build the table. On the
audit's synthetic 12-compound TCA-cycle input the Citrate cycle (hsa00020) p-value moved from
6.1e-19 (all-of-KEGG background, n=6709) to 9.4e-11 (320-ID assay-coverage background; re-run 2026-09-21) --
less significant with the correct, smaller background, as the theory predicts.

```bash
# reference_metabolome.txt: the assay-coverage background, one KEGG compound ID per line
Rscript scripts/local_ora.R kegg_ids.txt reference_metabolome.txt ora_local.csv   # columns: pathway, total, hits, p.value, fdr
```

### MetaboAnalystR API path (alternative; sends the compound list off-machine)

Use only if the remote call disclosed in Version Compatibility is acceptable. The script does its own compound mapping.

```bash
Rscript scripts/ora_api.R compounds.txt reference_metabolome.txt hsa name ora_api.csv   # columns include Raw p, FDR, Impact, Hits, Total
```

`scripts/ora_api.R` repeats the mapping steps (the API call needs the `mSet`), then `SetKEGG.PathLib` -> `SetMetabolomeFilter(mSet, TRUE)` -> `Setup.KEGGReferenceMetabolome()` -> `CalculateOraScore(mSet, 'rbc', 'hyperg')`.
`SetMetabolomeFilter(mSet, TRUE)` alone does NOT restrict the background: `Setup.KEGGReferenceMetabolome()` must be called first to load the reference file, or the filter silently has no effect (FALSE uses all of KEGG, the inflated default that manufactures false positives).
Checked on MetaboAnalystR 4.3.0: the server has been observed to reject the FILTERED request outright (`CalculateOraScore` returns 0; `current.msg == "Failed to connect to the API Server!"`), even with a correctly-matched reference file, while the unfiltered (FALSE) call to the same endpoint succeeds. The script then prints `current.msg` and exits 2; use Local-Only ORA above.
