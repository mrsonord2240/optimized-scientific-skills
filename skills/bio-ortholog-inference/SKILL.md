---
name: bio-ortholog-inference
category: Data Analysis
description: Pull pre-computed ortholog calls from public databases (OrthoDB, Ensembl Compara, OMA browser, eggNOG, PANTHER, KEGG Orthology, HomoloGene) via their REST APIs. Use when orthologs are already curated upstream, when the question is "what is the X ortholog of Y" rather than "how to infer orthology de novo", when batch-mapping gene IDs across species, or when comparing the resources for consensus calls. Encodes confidence-level semantics, 1:1 vs 1:many vs many:many, HomoloGene deprecation, and when to defect to de novo computation.
tool_type: python
primary_tool: requests
license: MIT
author: GPTomics
---

## Version Compatibility

Reference code needs requests 2.31+ and pandas 2.2+ (checked with 2.34.2 and 2.3.3). Resource versions: OrthoDB v12 API, Ensembl REST (release 116 live on 2026-10-02), OMA REST API, eggNOG 6.0+, PANTHER v18+.

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show requests pandas`
- API surface: confirm endpoint URLs and JSON schema match the current API docs

If endpoints return 404 or unexpected JSON, check release notes for the resource; schema migrations happen with each major version (Ensembl release is the biggest moving target).

# Ortholog Inference (Database Access)

**"What is the X ortholog of gene Y?"** -> Many ortholog resources have already done the inference at scale. Pulling their answers is faster and often more reliable than re-computing. This skill is the **database-access** view: how to query the major orthology resources programmatically, what their confidence semantics mean, and when their disagreements matter.

For **de novo orthology inference** (OrthoFinder, SonicParanoid, OMA standalone on local proteomes) see `comparative-genomics/ortholog-inference`.

Resources covered:
- **OrthoDB v12** -- broadest coverage (1700+ species), levels from species-specific to deep
- **Ensembl Compara** -- vertebrate-focused, tree-reconciled, one2one/one2many/many2many types with identity percentages
- **OMA browser** -- high precision, HOG (Hierarchical Orthologous Group) framework
- **eggNOG 6.0** -- pre-computed functional groups, deepest functional annotation
- **PANTHER** -- protein family + ortholog calls with experimentally validated curation
- **KEGG Orthology (KO)** -- pathway-centric orthologous functional units
- **HomoloGene** -- retired (see below)

No API keys are required for any of these (as of 2026), but rate limits apply.

## Client and examples

`scripts/ortholog_clients.py` is the single authoritative client (`sys.path.insert(0, 'scripts')`): `get_with_retry` (timeout; retries 429, 5xx and connection errors; used by every helper), `resolve_symbol`, `compara_orthologs`, `compara_one2one`, `batch_compara`, `orthodb_search`, `orthodb_groups`, `orthodb_orthologs`, `orthodb_group_has_gene`, `oma_orthologs`, `oma_hog_for_protein`, `oma_hog_members`, `ko_for_gene`, `genes_for_ko`, `ko_info`.

| Task | Example |
|---|---|
| Compara single gene + batch, symbol-rename check | `examples/compara_orthologs.py` |
| Compara vs OMA vs OrthoDB side by side for TP53 (per-resource error isolation, OrthoDB group verified) | `examples/cross_resource.py` |
| Gene -> KO -> members across species | `examples/kegg_orthology.py` |

```python
from ortholog_clients import compara_one2one, batch_compara

print(compara_one2one('Brca1', source='mouse', target='human'))   # dict or None; the Compara record has no confidence field
df = batch_compara(['BRCA1', 'TP53', 'MYC'], source='human', target='mouse')
print(df[df['type'] == 'ortholog_one2one'][['symbol', 'target_id', 'taxonomy_level', 'target_pid']])
```

## Decision matrix: which resource for which question?

| Question | Resource | Why |
|---|---|---|
| Ortholog of human gene X in mouse | Ensembl Compara | Best-curated for vertebrates; ortholog type and identity per call |
| Ortholog of gene X across all 1700+ species | OrthoDB | Broadest taxonomic coverage |
| Single-copy orthologs for phylogenomics | OrthoDB at species-tree level | Pre-computed; large taxonomic groups |
| Functional annotation transfer | eggNOG-mapper | OG-based functional categories |
| Pathway-centric orthology | KEGG Orthology (KO) | KO IDs link directly to pathway maps |
| Curated function-aware orthologs | PANTHER | Smaller scope; manually curated; experiment-supported |
| Compare resource consensus | All of them + intersect | Disagreement is itself a signal |
| Plant orthology | Ensembl Compara (plants division) | Better than the vertebrate set for plants |
| Bacterial orthology | OrthoDB or eggNOG bactNOG | Ensembl Bacteria has limited Compara coverage |
| Custom proteomes not in any database | **De novo computation** | See `comparative-genomics/ortholog-inference` |

## Per-resource API reference

### OrthoDB v12

Base `https://data.orthodb.org/v12/`, GET, JSON:
- `/search?query=<text>&species=<NCBI_taxid>` -- **full-text** search over group descriptions, not a gene-symbol lookup. Returns `{'data': [og_id, ...], 'bigdata': [{id, name, level_name, gene_count, ...}]}` in relevance order; `TP53` ranks phosphoglycerate mutase (TIGAR) and TP53-target groups first and `BRCA1` ranks Ubiquitin first, while `tumor protein p53` reaches the p53 group 4289813at2759. Verify a group (name, or `orthodb_group_has_gene`) before use.
- `/orthologs?id=<og_id>&species=<taxid>` -- `data` is a list of groups, each with `genes[].gene_id.id` (symbol); `data` is `null` when the group has no member at that level
- `/group?id=<og_id>` -- full group info (name, level, InterPro domains)
- `/tab?id=<og_id>` -- tab-separated member table (the parameter is `id`; `query=` returns 404)

Levels are NCBI taxonomy IDs (9606 human, 40674 Mammalia, 7742 Vertebrata). Bulk downloads are available.

### Ensembl Compara (via Ensembl REST)

Base `https://rest.ensembl.org/`, `Accept: application/json`, 15 req/sec and 55,000 req/hour, `Retry-After` on 429.
- `/homology/symbol/<species>/<symbol>` or `/homology/id/<species>/<ensembl_gene_id>`
- `/lookup/symbol/<species>/<symbol>` -- resolve symbol to Ensembl ID first
- `?type=orthologues` drops paralogs; `?target_species=<species>` restricts to one target
- Types: `ortholog_one2one`, `ortholog_one2many`, `ortholog_many2many`, `within_species_paralog`

For broader Ensembl workflows see `ensembl-rest`.

### OMA REST API

Base `https://omabrowser.org/api/`, JSON, no published rate limit (be polite):
- `/protein/<id>/orthologs/` -- orthologs of a protein (UniProt or OMA ID); add `?rel_type=1:1` to filter server-side (the most reliable call; unfiltered calls 502 intermittently). Items carry `omaid`, `canonicalid`, `species.code`, `species.taxon_id`. A 200 with `[]` means OMA has no orthologs (BRCA1 P38398), not an outage; TP53 `P04637` returns `P53_MOUSE`
- `/protein/<id>/` -- protein record including `oma_hog_id`
- `/hog/<hog_id>/` -- Hierarchical Orthologous Group info
- `/genome/<species_code>/` -- genomes; species codes are 5-letter (HUMAN, MOUSE)

### eggNOG 6.0

The eggNOG web API was not usable on 2026-10-03: `eggnogdb.org/api` returned 403 to a scripted client and `eggnog6.embl.de` failed TLS verification; treat it as unverified. For batch protein-set annotation use **eggNOG-mapper** locally (Cantalapiedra et al. 2021 *Mol Biol Evol* 38:5825); for ad hoc lookup search the eggNOG web interface for an orthogroup ID, then download the member set.

### KEGG Orthology (KO)

Base `https://rest.kegg.jp/`; returns plain-text TSV, **not** JSON (an HTML reply means a wrong endpoint). `ko_for_gene('hsa', '7157')` maps human TP53 to its KO; `genes_for_ko` lists all member genes across KEGG species. KEGG license: free for academic web/API use, paid for commercial use.

### PANTHER

Base `https://pantherdb.org/services/oai/pantherdb/` (the `http://` form redirects to https). Smaller curated scope than OrthoDB. The client does not wrap PANTHER. Verified route: `GET ortholog/matchortho?geneInputList=P04637&organism=9606&targetOrganism=10090&orthologType=LDO` returns the target gene (mouse Tp53, UniProtKB P02340) and an ortholog type code (`LDO`, least-diverged ortholog); the response carries no evidence-code field (checked 2026-10-03, PANTHER v19).

### HomoloGene (retired)

Frozen since 2014, and NCBI E-utilities no longer serve it: `efetch.fcgi?db=homologene` returned "Database ... not supported" on 2026-10-02. Treat as legacy, do not build live workflows on it, and verify any historical HomoloGene call against Compara or OrthoDB.

## Confidence-level semantics

Each resource defines "confidence" differently; do not average or compare across resources. Use within-resource cutoffs and intersect call sets for cross-resource comparison.

| Resource | Confidence field | Semantics |
|---|---|---|
| Ensembl Compara | none in live records (verified 2026-10-03) | Use `type`, `taxonomy_level`, `perc_id` of source and target; `compara_orthologs` still returns `confidence` as None |
| OrthoDB | `evolutionary_rate` (not a confidence) | Inverse proxy; lower = more conserved |
| OMA | Internal QC; not exposed per call | All calls passed a precision filter |
| eggNOG | Tax-level coverage | Member counts per taxonomic level |
| PANTHER | ortholog type (e.g. `LDO`) | Least-diverged vs other ortholog; no per-call evidence field in the response |

## The orthology conjecture (and why resources disagree)

Orthologs are more likely than paralogs to share function (Tatusov 1997), but only weakly (Studer & Robinson-Rechavi 2009 *Trends Genet* 25:210; Altenhoff et al. 2012 *PLoS Comput Biol* 8:e1002514); sub- and neo-functionalization let a paralog become the functional equivalent. Resources disagree because their algorithms weigh evidence differently:
- **OMA** is strict (RBH + verification + HOG inference): higher precision, lower recall.
- **Ensembl Compara** is tree-reconciled: best for vertebrates.
- **OrthoDB** uses broad hierarchical clustering: wider coverage, more ambiguous calls.
- **eggNOG** uses pre-computed orthogroups at fixed taxonomic levels: fast but coarser.

For high-stakes calls (publication, drug-target choice) **intersect at least two resources** after normalizing IDs to a common namespace (UniProt or NCBI Gene; see `uniprot-access`) and inspect the disagreements.

## Scale

REST loops are fine for hundreds of genes (`batch_compara`, 0.07 s sleep); above ~5,000 genes use Ensembl BioMart bulk export (`biomart-queries`).

## Failure modes

- **Resource disagreement on 1:1** -- Compara 1:1, OrthoDB 1:many, OMA no call. Define the authoritative resource per project or take the intersection; document the choice.
- **Stale snapshot** -- a 2-year-old OrthoDB download or HomoloGene misses current orthologs and can cite defunct IDs. Pin a dated release; refresh annually.
- **Symbol ambiguity** -- `MARCH1` was renamed `MARCHF1` in 2020; lookups return empty or the wrong gene. Resolve to the Ensembl/HGNC ID first.
- **No Compara `confidence`** -- do not filter on it; filter on `type` (`ortholog_one2one`) and identity.
- **OrthoDB first hit is not the gene** -- `/search` is full text; verify the group before reading its members.
- **Rate-limit cascade** -- 5000 Ensembl calls without sleep give 429 stalls; sleep and honor Retry-After.
- **Custom proteome not in any DB** -- no calls are returned; run OrthoFinder/SonicParanoid (see `comparative-genomics/ortholog-inference`).
- **KEGG license confusion** -- commercial use needs a paid license; eggNOG and OrthoDB are more permissive.

## Common errors

| Error / symptom | Cause | Solution |
|---|---|---|
| `HTTPError 429` (Ensembl) | Rate limit | Sleep per Retry-After; cap at 15 req/sec |
| Empty homologies | Symbol misspelled or stale | Resolve to Ensembl ID first |
| Compara `confidence` always None | Field absent from current Ensembl records | Use type, taxonomy level and identity |
| OMA `404` | Wrong namespace (UniProt where OMA ID needed) | Use the protein lookup endpoint to resolve first |
| OMA `502` | Transient gateway failure | The client retries with backoff; prefer `rel_type=1:1`; if it persists the helper raises, report the leg as unavailable |
| KEGG returns HTML | Wrong endpoint | Use `rest.kegg.jp`; output is TSV not JSON |
| Resource disagreement | Different algorithms / coverage | Intersect; document choice |

## References

- Tatusov RL, Koonin EV, Lipman DJ. (1997) A genomic perspective on protein families. *Science* 278:631-637.
- Altenhoff AM, Studer RA, Robinson-Rechavi M, Dessimoz C. (2012) Resolving the ortholog conjecture: orthologs tend to be weakly, but significantly, more similar in function than paralogs. *PLoS Comput Biol* 8:e1002514.
- Studer RA, Robinson-Rechavi M. (2009) How confident can we be that orthologs are similar, but paralogs differ? *Trends Genet* 25:210-216.
- Kuznetsov D, Tegenfeldt F, Manni M, Seppey M, Berkeley M, Kriventseva EV, Zdobnov EM. (2023) OrthoDB v11: annotation of orthologs in the widest sampling of organismal diversity. *Nucleic Acids Res* 51:D445-D451.
- Herrero J, Muffato M, Beal K, et al. (2016) Ensembl comparative genomics resources. *Database* 2016:baw053.
- Altenhoff AM, Vesztrocy AW, Bernard C, et al. (2024) OMA orthology in 2024. *Nucleic Acids Res* 52:D513-D521.
- Hernandez-Plaza A, Szklarczyk D, Botas J, et al. (2023) eggNOG 6.0: enabling comparative genomics across 12,535 organisms. *Nucleic Acids Res* 51:D389-D394.
- Cantalapiedra CP, Hernandez-Plaza A, Letunic I, Bork P, Huerta-Cepas J. (2021) eggNOG-mapper v2: functional annotation, orthology assignments, and domain prediction at the metagenomic scale. *Mol Biol Evol* 38:5825-5829.
- Thomas PD, Ebert D, Muruganujan A, Mushayahama T, Albou LP, Mi H. (2022) PANTHER: making genome-scale phylogenetics accessible to all. *Protein Sci* 31:8-22.

## Related Skills

- comparative-genomics/ortholog-inference - De novo orthology computation (OrthoFinder, OMA standalone, SonicParanoid)
- ensembl-rest - Broader Ensembl REST workflows beyond Compara
- biomart-queries - Bulk ortholog table export via Ensembl BioMart
- uniprot-access - Resolve UniProt accessions used by OMA orthology queries
- pathway-analysis/kegg-pathways - KEGG Orthology and pathway mapping
