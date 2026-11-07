# Explore WikiPathways

```bash
Rscript examples/wikipathways_explore.R
```

This illustrative script contacts the live API and archive, lists organisms and pathways, searches text, retrieves metadata and cross-references, and demonstrates a dated GMT. It uses hard-coded human example values and an archive-derived universe for its demonstration: adapt those values and supply the actual tested universe for an analysis. Search and metadata calls may be extracted from the example; there is no general parameterized query CLI. `listOrganisms()` supplies exact organism names; cross-reference system codes include L (Entrez), H (HGNC), and En (Ensembl). Inspect pathway metadata and last-edit date before relying on a hit.

The retained `scripts/wikipathways_api_smoke.R` is a source-provided live API regression example, not a general query interface; it was not executed during normalization.

**Done when:** the selected organism, identifiers, metadata and database release or retrieval date are recorded.
