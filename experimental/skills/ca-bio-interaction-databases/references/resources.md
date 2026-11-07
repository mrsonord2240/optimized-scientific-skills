# Resources



| Question | Best resource | Why |
|---|---|---|
| "Build a network around 10 genes" | STRING (medium confidence ~400) | Comprehensive; channels combinable; good viz integration |
| "Only physically interacting proteins" | IntAct or BioGRID physical (or STRING `network_type=physical`) | Curated physical interactions; PSI-MI standard |
| "Signed/directed signaling (phospho, ubiq, etc.)" | **SIGNOR** | Only major DB with mechanism types and direction |
| "Functional enrichment based on co-mentioned genes" | STRING functional (default) | Includes textmining channel |
| "Genetic interactions (synthetic lethality)" | BioGRID genetic | Largest curated genetic interaction set |
| "High-throughput Y2H interactome" | **HuRI** | Reference yeast-2-hybrid map of human |
| "Mass-spec-derived protein complexes" | **HuMAP v2** or BioPlex | AP-MS complex maps |
| "Curated pathways with interactions" | Reactome | Pathway-organized; gold standard for signaling |
| "Meta-database aggregating 100+ sources" | **OmniPath** | The modern "one-stop"; pre-aggregated |
| "Cross-species or non-human" | STRING | Species coverage broadest |
| "Bacterial interactome" | STRING bacterial | Limited curated alternatives |
| "Phosphorylation site-specific" | PhosphoSitePlus (commercial license) or SIGNOR | PSP has best PTM coverage but requires license |

## STRING (v12.x, channels, confidence)

STRING (Szklarczyk et al. 2023 *Nucleic Acids Res* 51:D638) aggregates evidence into a combined confidence score over seven **channels**:

| Channel | TSV column | What it captures |
|---|---|---|
| `experiments` | `escore` | Direct experimental evidence (BioGRID, IntAct, etc.) |
| `database` | `dscore` | Curated database (Reactome, KEGG) |
| `textmining` | `tscore` | Co-mention in PubMed |
| `coexpression` | `ascore` | Co-expression across conditions |
| `neighborhood` | `nscore` | Genomic neighborhood (prokaryotes mainly) |
| `fusion` | `fscore` | Gene fusion across species |
| `cooccurrence` | `pscore` | Phylogenetic profile co-occurrence |

`required_score` is the combined score on a 0-1000 scale; TSV output reports 0-1 scores. `string_network(..., network_type='physical')` returns STRING's physical-complex network; the `experiments` channel (`escore`) mixes physical and functional experimental evidence, so `escore > 0.4` is not a physical filter (for 10 DNA-damage genes at 700, 5 of 26 such edges were absent from the physical network).

| Threshold | Tier | Use when |
|---|---|---|
| 150 | Low | Exploration; includes weak textmining |
| 400 | Medium | Default; balanced sensitivity/specificity |
| 700 | High | Publication-quality networks |
| 900 | Highest | Experimentally validated core only |

`caller_identity=<app-name>` lets STRING attribute and contact heavy automated callers; non-compliant clients get throttled first. STRING wants one request at a time per `caller_identity`: sleep 1-2 s between calls, and use bulk downloads for very large batches.

## BioGRID (physical + genetic, curated + HT)

BioGRID (Oughtred et al. 2021 *Protein Sci* 30:187) covers physical and genetic interactions across the broadest organism set; free **API key** required. Key fields: `EXPERIMENTAL_SYSTEM` (e.g. "Two-hybrid", "Affinity Capture-MS", "Synthetic Growth Defect"), `THROUGHPUT` ("Low Throughput" vs "High Throughput" -- the most important quality flag), `EVIDENCE_TYPE` (physical vs genetic). For high-confidence physical interactions keep `THROUGHPUT == 'Low Throughput'` and `EXPERIMENTAL_SYSTEM in PHYSICAL_LT_SYSTEMS` (`biogrid_lt_physical`). Use the NCBI taxon ID for `taxId`. A missing or rejected key returns HTTP 401 (raised by the client without the URL); an empty result means no matching rows, usually a wrong taxon or gene name. `biogrid_lt_physical` returns one row per experiment, not per pair (TP53: 3081 rows for 831 pairs); deduplicate on the sorted gene pair before counting interactions.

## IntAct (PSI-MI curated physical)

IntAct (Del Toro et al. 2022 *Nucleic Acids Res* 50:D648) is the IMEx consortium reference for curated physical interactions in PSI-MI format. MINT was folded in ~2014; MINT URLs redirect to IntAct. The direct REST API has changed multiple times; the most stable access is via OmniPath (which wraps IntAct) or the PSICQUIC web services.

## SIGNOR (signed signaling)

SIGNOR (Lo Surdo et al. 2023 *Nucleic Acids Res* 51:D631) is **the only major curated database with signed, directed, mechanism-typed interactions**. Each edge has `direction` (A->B), `effect` (`up-regulates`, `down-regulates`, `unknown`) and `mechanism` (`phosphorylation`, `dephosphorylation`, `ubiquitination`, `binding`, ...). Essential for signaling-pathway analysis and dynamic modeling; coverage is smaller than STRING/BioGRID but quality is high.

`signor_for_gene(symbol)` resolves the symbol to a reviewed human UniProt accession (UniProt REST; `uniprot_accession`) and calls `getData.php?organism=9606&id=<accession>`. The response is headerless with 29 tab-separated columns (verified 2026-10-03): ENTITYA col 0, ENTITYB col 4, EFFECT col 8, MECHANISM col 9, residue col 10, PMID col 21, SIGNOR id col 26, score col 27. `organism=human` returns "No result found." and `entity=<symbol>` is ignored. Rows are per record (site/paper), cover the protein as either end, and `effect` is free text such as `up-regulates quantity by stabilization`. TP53 returned 333 rows. 'No result found.' yields an empty frame plus a warning; an unparseable non-empty answer raises ValueError.

## Reactome (curated pathways)

Reactome (Milacic et al. 2024 *Nucleic Acids Res* 52:D672) is the gold standard for human pathway curation with full interaction reactions; species coverage outside human is limited. REST: `https://reactome.org/ContentService/`.

## HuRI / HuMAP (human-specific interactomes)

- **HuRI** (Luck et al. 2020 *Nature* 580:402): yeast-2-hybrid interactome of human; ~53K binary interactions; biased toward binary high-confidence. Portal: `http://www.interactome-atlas.org/`.
- **HuMAP v2** (Drew et al. 2021 *Mol Syst Biol* 17:e10016): AP-MS-derived complex map integrating >15,000 proteomic experiments.

Both are available for download.

## OmniPath (meta-database)

OmniPath (Türei et al. 2021 *Mol Sys Biol* 17:e9923) aggregates 100+ sources with provenance; the modern default for "everything anyone has said about A-B".

| Endpoint | Content |
|---|---|
| `/interactions` | Signaling interactions (directed); most useful for cross-DB consensus |
| `/enzsub` | Enzyme-substrate (kinase-substrate, etc.) |
| `/complexes` | Protein complexes |
| `/annotations` | Functional annotations |
| `/intercell` | Intercellular communication |

Each interaction carries `sources` and `references` (PMIDs). Parameters: `types=` (post_translational, transcriptional, mirna_target, lncrna_target), `license=academic|commercial`. The server accepts `license` but does not filter on it (TP53 gave identical rows, including PhosphoSite and HPRD, which OmniPath lists as academic-purpose), so the client screens client-side for `commercial` using each resource's license purpose from `/resources`: sources not marked commercial (academic, non-profit, composite, unknown) are dropped and rows with no permitted source are removed. This is a source-level screen only; per-resource terms (NC, ND, share-alike) and server-computed fields still need review. `omnipath_interactions` returns each gene's whole neighbourhood, not only edges among the query genes. R users: `OmnipathR`.



- `aggregate_networks(genes, biogrid_key=None)` builds an undirected union of edges among the query genes (OmniPath and BioGRID neighbour edges are restricted to the query set so scopes match) with a `sources` set and `string_score` (0-1, None when STRING lacks the edge) per edge (self-interactions, i.e. homodimers, are excluded from the aggregate, so a gene whose only evidence is a self-interaction is not a node; those `gene_a == gene_b` rows remain in `biogrid_lt_physical` output); `multi_source_edges(g)` returns edges reported by two or more resources (overlap, not calibrated confidence); `summary(g)` gives nodes/edges/density/components/mean degree.
- SIGNOR and OmniPath edges are directional and signed: use `nx.DiGraph` and keep `effect` and `mechanism`. STRING edges are undirected and require both directions if represented in a directed graph.
- Resolve gene symbols to UniProt/HGNC IDs before querying when possible; symbols drift (MARCH1 -> MARCHF1 in 2020).

## Failure modes

- **Functional treated as physical** -- STRING's combined score mixes textmining and coexpression. `escore` still mixes physical and functional experiments; use `network_type=physical`, or BioGRID/IntAct, for physical claims.
- **Wrong confidence tier** -- default 400 in a publication network includes weak textmining hits; use 700 (900 for the validated core).
- **HT trusted as LT** -- Y2H/AP-MS screens have higher false-positive rates; filter `THROUGHPUT = 'Low Throughput'` or use IntAct curated scores.
- **Symbol drift** -- renamed symbols give empty or wrong matches; resolve to IDs first.
- **STRING version drift** -- code pinned to `version-11-5` silently returns stale data; pin the current version or use `string-db.org`.
- **License surprise** -- ConsensusPathDB is academic-only; audit each source; OmniPath `license=commercial` is screened client-side and still needs review.
- **Direction collapsed** -- treating SIGNOR/OmniPath edges as undirected loses topology; use DiGraph.
- **STRING rate limiting** -- looping `/network` over many gene sets triggers throttling; sleep 1-2 s or use bulk downloads.

## Common errors

| Error / symptom | Cause | Solution |
|---|---|---|
| Stale STRING data | Old pinned `version-11-5` host | Use the current version or unversioned host |
| BioGRID HTTP 401 | Missing or rejected API key | Get a key; never print it |
| BioGRID empty result | Wrong taxId or gene name | Use NCBI taxon ID and official symbol |
| Symbol mismatch | HGNC renaming | Resolve via UniProt/Ensembl ID |
| HT interactions inflate network | No throughput filter | Filter `THROUGHPUT = 'Low Throughput'` |
| Functional vs physical confusion | Mixed STRING channels | `network_type=physical` |
| Directional edges collapsed | Used Graph for directed source | Use DiGraph for SIGNOR/OmniPath |
| License violation in commercial pipeline | ConsensusPathDB or PhosphoSitePlus | Switch to permissive sources |
| OmniPath returns nothing | `license=commercial` screen removed every source | Use `license='academic'` for academic use |
| SIGNOR empty frame and warning | Symbol has no SIGNOR record, or wrong accession | Check `uniprot_accession(symbol)` |
