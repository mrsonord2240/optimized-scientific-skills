# Endpoints



Ensembl Genomes was folded into the main REST service: `https://rest.ensembl.org` serves every division (`/info/divisions` lists EnsemblVertebrates, EnsemblPlants, EnsemblFungi, EnsemblMetazoa, EnsemblProtists, EnsemblBacteria). `rest.ensemblgenomes.org` no longer resolves (checked 2026-10-02).

| Need | Call |
|---|---|
| Enumerate species of a division | `/info/species?division=EnsemblPlants` |
| Plants / fungi / metazoa / protists / bacteria | same host; use the species name from `/info/species`, e.g. `/lookup/symbol/arabidopsis_thaliana/NAC001` |
| Bacteria | limited Compara coverage; most bacterial work belongs in NCBI |

## Version pinning

| URL | Behavior |
|---|---|
| `https://rest.ensembl.org` | Current release (rolling) |
| `https://e110.rest.ensembl.org` | Pinned to release 110 (redirects to the monthly archive host) |
| `https://e111.rest.ensembl.org` | Pinned to release 111 |
| `https://grch37.rest.ensembl.org` | Pinned to GRCh37 (legacy assembly) |

For any published analysis, pin the release. Ensembl releases change gene model versions, exon coordinates, and transcript annotations -- re-running a pipeline a year later against the live endpoint may produce different results. Re-runs against an archive host return the same Gene ID even after the live release moves on. Record the pinned release in pipeline metadata. Old archives are retired roughly five years after release (DNS failure or 503); see https://www.ensembl.org/info/website/archives/index.html.

## Major endpoint groups

| Group | Example | Purpose |
|---|---|---|
| Lookup | `/lookup/symbol/human/BRCA1` | Resolve symbol or ID to stable record |
| Sequence | `/sequence/id/{id}` | DNA/protein sequence for ID |
| Cross References | `/xrefs/symbol/human/BRCA1` | Cross-refs to other DBs |
| Homology / Compara | `/homology/symbol/human/BRCA1` | Orthologs and paralogs |
| Gene Tree | `/genetree/id/{tree_id}` | Compara gene tree |
| VEP | `/vep/human/region/{region}/{allele}` | Variant effect prediction |
| Overlap | `/overlap/id/{id}` or `/overlap/region/{species}/{region}` | Genes/regulatory in interval |
| Regulatory | `/overlap/region/{species}/{region}?feature=regulatory` | Regulatory features in an interval (no verified per-ID route; `/regulatory/...` ID paths returned 404) |
| Variant | `/variation/{species}/{id}` | dbSNP / 1000G / ClinVar via Ensembl |
| LD | `/ld/{species}/pairwise/{var1}/{var2}` | LD between variants |
| GA4GH | `/ga4gh/...` | GA4GH-compliant subset |

Full reference: https://rest.ensembl.org (interactive). Species names are lowercase Ensembl names or common names (`human`, `homo_sapiens`); `Homo_sapiens` was also accepted when checked 2026-10-02; prefer lowercase.

## The symbol-vs-ID stability problem

Gene symbols are unstable (MARCH1 -> MARCHF1 in 2020 due to Excel autocorrect; SEPT* family also renamed). Ensembl Gene IDs (ENSG...) are stable across releases when the gene model is preserved.

1. Resolve symbol -> Ensembl ID once at pipeline start: `/lookup/symbol/{species}/{symbol}`.
2. Persist the Ensembl ID.
3. Run downstream queries by ID, not symbol.

Symbol endpoints suit interactive use; ID endpoints suit reproducible pipelines. A renamed or unknown symbol returns HTTP 400 (`No valid lookup found`) or, rarely, the wrong gene.

## Rate limits

| Limit | Value |
|---|---|
| Burst | 15 req/sec |
| Hourly | 55,000 req/hour |
| Concurrent | Not enforced; courtesy 1-2 |

Sleep 0.07 s between calls (`SLEEP`), honor `Retry-After` on 429 (`get_with_retry` does). For >5,000 queries switch to BioMart bulk export (see `bio-biomart-queries`). Some endpoints accept batches in a POST body (e.g. `/lookup/id` with `{"ids": [...]}`); use that for <500 IDs to cut request count.

## VEP (Variant Effect Predictor)

VEP via REST is for ad hoc annotation (<1K variants). For batch work (>1000 variants) run VEP locally (`variant-calling/variant-annotation` skill); 100K REST calls take days and trigger 429 cascades.

- `/vep/{species}/region/{region}/{allele}` -- by coordinate, `chr:start-end:strand`
- `/vep/{species}/id/{variant_id}` -- by dbSNP / Ensembl variant ID
- `/vep/{species}/hgvs/{hgvs_notation}` -- by HGVS notation

Returns consequence terms, SIFT/PolyPhen, gnomAD frequencies and ClinVar significance where available. `rest.ensembl.org` defaults to GRCh38; for GRCh37 coordinates use `https://grch37.rest.ensembl.org`.

## Compara homology

Detailed orthology resource comparison lives in `bio-ortholog-inference`. Endpoints:

- `/homology/symbol/{species}/{symbol}` or `/homology/id/{species}/{ensembl_id}` (species segment required) -- all orthologs across Ensembl species
- `?target_species=` restricts to one target; `?type=orthologues` or `paralogues` filters
- Types seen: `ortholog_one2one`, `ortholog_one2many`, `ortholog_many2many`, `within_species_paralog`; `confidence` is often absent (treat missing as unknown). Paralog lists can be empty (BRCA1 returned none)

## Failure modes

- **Symbol pipeline breaks on HGNC rename** -- HTTP 400 or wrong gene (e.g. `MARCH1` after 2020). Resolve to Ensembl Gene ID once; use IDs downstream.
- **No version pinning, results drift** -- live endpoint follows the current release, so coordinates, exon counts, sometimes Gene IDs change. Pin an archive endpoint.
- **Rate-limit cascade** -- a loop of 5,000 calls without sleep stalls on 429. Sleep 0.07 s, honor Retry-After, use BioMart above 5K.
- **VEP for bulk variants** -- infeasible over REST; run VEP locally.
- **Wrong species name** -- 404 or empty result; enumerate valid names with `/info/species`.
- **Non-vertebrate species not found** -- use the exact species name from `/info/species?division=...`; do not use the retired `rest.ensemblgenomes.org` host.
- **Archive TLS / connectivity** -- very old archives (e80, e90) are decommissioned (DNS failure, 503, or a 200 HTML page that the client rejects as non-JSON); a recent archive host can also be down temporarily (e116 returned 503 on 2026-10-03 while the live host worked), so retry later or use another release.

## Common errors

| Error / symptom | Cause | Solution |
|---|---|---|
| `400` on symbol lookup (`No valid lookup found`) | HGNC rename, unknown symbol or wrong species | Resolve to Ensembl ID; check species code |
| HTTP 429 | Rate limit | Sleep per Retry-After; cap to 15 req/sec |
| Different results 6 months later | No version pinning | Use archive endpoint (`eXX.rest.ensembl.org`) |
| `404` on non-vertebrate | Wrong species name | Use `/info/species?division=...` |
| VEP infeasible bulk | REST is per-variant | Local VEP for bulk |
| Old archive 503 / HTML page | Decommissioned, or a recent host temporarily down | Retry later; use another release |
| `400` on `/sequence/id/ENSG...?type=protein` | Gene ID with a non-genomic type | Use an ENST/ENSP ID or `multiple_sequences=True` |
| `400` on VEP region | Coordinate from the wrong assembly | GRCh37 coordinates go to `grch37.rest.ensembl.org` |
