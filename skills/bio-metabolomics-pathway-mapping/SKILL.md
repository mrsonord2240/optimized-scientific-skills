---
name: bio-metabolomics-pathway-mapping
category: Data Analysis
description: Maps metabolomics results to biological pathways via over-representation (ORA), metabolite-set enrichment (MSEA/QEA), mummichog/PSEA on raw m/z peaks, and network-diffusion enrichment (FELLA), with correct background-set construction and honest interpretive ceilings. Use when interpreting differential metabolites or an untargeted LC-MS feature table in pathway context, choosing ORA vs MSEA vs mummichog vs topology, or setting the reference/background set. For annotation confidence levels feeding ORA see metabolomics/metabolite-annotation; for gene-set concepts see pathway-analysis/go-enrichment and pathway-analysis/gsea; for joint gene+metabolite pathways see multi-omics-integration/mofa-integration.
tool_type: r
primary_tool: MetaboAnalystR
license: MIT
author: GPTomics
---

## Version Compatibility

Reference examples tested with: MetaboAnalystR 4.0+, FELLA 1.22+

Install: MetaboAnalystR from GitHub (https://github.com/xia-lab/MetaboAnalystR); `BiocManager::install(c("FELLA", "KEGGREST"))`.
`PerformPSEA` also needs `install.packages(c("fitdistrplus", "RJSONIO"))` -- declared only under
MetaboAnalystR's Suggests, so a `dependencies=FALSE` GitHub install omits them and PSEA throws
`there is no package called 'RJSONIO'` (checked on MetaboAnalystR 4.3.0).

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

The single most important input fact: whether the metabolites are confidently identified (KEGG/HMDB IDs) determines which method is even possible. An untargeted LC-MS feature table with no IDs cannot run ORA; it requires mummichog/PSEA. Verify the input type before choosing a tool.

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

**Required local-session setup (checked on MetaboAnalystR 4.3.0).** Outside the MetaboAnalyst
Shiny web app, MetaboAnalystR's own error handler (`AddErrMsg`) does `current.msg <<- c(current.msg, msg)`
against globals that only the web app initializes; in a plain `library(MetaboAnalystR)` session this
crashes with `object 'current.msg' not found` on *any* local validation failure (bad background, too
few features, bad mapping) -- hiding the real diagnostic message behind an unrelated R error. Before
any MetaboAnalystR call in this Skill, run:
```r
current.msg <- character(0); err.vec <- character(0)
```
If a call then returns `0` instead of an `mSet` object, print `current.msg` for the actual reason.

**Undisclosed remote call.** `CalculateOraScore()` and `CalculateQeaScore()` against a KEGG pathway
library do not compute locally, even when run outside the website: `my.ora.kegg()` serializes the
mapped compound list and POSTs it (`httr::POST`, multipart, serialized `mSet`) to
`https://www.xialab.ca/api/pathwayora` (or `/pathwayqea`) for every call, filtered or not -- confirmed
by inspecting `my.ora.kegg`/`.do.api.call`'s source. A user running "locally installed" MetaboAnalystR
on unpublished compound names/IDs should know they are leaving the machine before running this. If
that is not acceptable, use the **Local-Only ORA** function (`references/ora-local-and-api.md`) instead, which never leaves the
machine (it only downloads the public, organism-agnostic KEGG pathway-to-compound reference table --
no user data is sent).

**Reference-library downloads (not user data).** `SetKEGG.PathLib()`, `CrossReferencing()`,
`Setup.KEGGReferenceMetabolome()` and the mummichog/`PerformPSEA` path load their libraries through
MetaboAnalystR's internal `.get.my.lib(filenm, sub.dir)`. Outside the web app it runs
`download.file("https://www.metaboanalyst.ca/resources/libs/...")` into the *working directory*
whenever the cached copy is missing or older than 30 days (printing `"Loaded files from MetaboAnalyst
web-server."`). Only generic compound/pathway/adduct libraries are fetched -- the function never
receives the compound list or `mSet` (checked by reading its source on 4.3.0) -- but the first run
needs internet, results can shift when a stale cache refreshes, and an offline machine must pre-seed
those files.

# Metabolomics Pathway Mapping

**"Map my metabolites to pathways"** -> Test whether a metabolite set or an m/z feature table is statistically enriched for biochemical pathways, given an explicit background.
- Identified compound list -> ORA / MSEA: `CalculateOraScore()` (MetaboAnalystR)
- Raw m/z peak table (no IDs) -> mummichog / PSEA: `PerformPSEA()` (MetaboAnalystR)
- Mechanism (which enzymes/reactions link the hits) -> network diffusion: `runDiffusion()` (FELLA)

## The Single Most Important Modern Insight -- Pathway Enrichment Launders Annotation Uncertainty Into Confident Biology

Enrichment is the one workflow step where uncertainty is structurally destroyed: compounds enter as names with no error bars, and the hypergeometric/permutation machinery cannot represent "this is a 40%-confident guess." A pile of MSI-level-3 tentative annotations emerges as a p-value with three decimals. Wieder 2021 simulated this directly: even a 4% misidentification rate manufactured both false-positive and false-negative pathways across five real datasets, and real untargeted annotation is far worse than 4%. The errors do not average out, because a single wrong hub-adjacent compound (alanine, glutamate, a TCA intermediate) can flip a pathway by itself. No untargeted pathway claim can be stronger than its annotation layer. The honest ceiling is "features consistent with perturbation of pathway X co-varied with phenotype, conditional on the chosen annotations, background, database boundary, and ionization settings" -- never "pathway X is upregulated."

## Two Starting Points -> Method -> Tool

The field's most common category error is conflating identified-compound enrichment with raw-feature activity prediction. They are disjoint entry points.

| Input | Goal / situation | Method | Tool | Key constraint |
|---|---|---|---|---|
| Identified compounds + cutoff | Discrete "significant" hit list | ORA (hypergeometric) | MetaboAnalystR `CalculateOraScore` | Background = compounds the assay could detect, NOT all of KEGG; code in `references/ora-local-and-api.md` |
| Identified compounds + ranked stat | No natural cutoff; keep magnitude | MSEA / QEA (rank-aware) | MetaboAnalystR `CalculateQeaScore` | Needs a meaningful, complete ranking |
| Raw m/z + RT + per-feature stat, NO IDs | Predict pathway activity, bypass ID | mummichog / GSEA-PSEA | MetaboAnalystR `PerformPSEA` | Background = the FULL feature table (R_all); declare ionization mode; code in `references/mummichog-psea.md` |
| Identified compounds | Mechanism: which enzymes/reactions link hits | Network diffusion | FELLA `runDiffusion` | KEGG IDs only; check `getExcluded()` for unmapped; code in `references/fella-diffusion.md` |
| Identified compounds | Database coverage is the bottleneck | Chemical-structure clustering | ChemRICH (background-independent) | Sidesteps pathway dark matter |
| Any | Secondary lens only | Topology / "impact" | MetaboAnalystR (MetPA) | Hub artifact; never sole evidence |

Mummichog exists because identification is the rate-limiter: only ~2-10% of untargeted features are ever confidently identified. It predicts network activity directly from the feature table, then the network context retro-prioritizes which annotation was probably right (Li 2013). Its existence is an admission of the annotation bottleneck, not a triumph -- use it knowing systems-level inference is bought with per-metabolite certainty.

## ORA vs MSEA vs Topology vs Mummichog

| Axis | ORA (hypergeometric) | MSEA / GSEA-PSEA | Topology / "Impact" | Mummichog / PSEA |
|---|---|---|---|---|
| Input | Identified list + cutoff | Identified ranked list | Identified list in pathway graphs | Raw m/z + RT + stat, no IDs |
| Null question | More hits than chance? | Set systematically high/low in ranking? | Hits at central (high-betweenness) nodes? | Do mass-matched candidates cluster in pathways beyond a random feature list? |
| Uses magnitude? | No (cutoff discards it) | Yes | Indirectly (enrichment x centrality) | No (cutoff defines the query) |
| Null source | Assay-coverage background | The ranked universe | Curated graph structure | Permutation from the FULL feature table (R_all) |
| Headline failure | Wrong/implicit background | Needs a complete ranking | Hub overemphasis (alanine ~95% case) | Significant-features-only as background |
| Output | Measured enrichment | Measured enrichment | Graph property, not the experiment | PREDICTED activity, not identities |

## ORA on an Identified Compound List

**Goal:** Test whether a list of confidently identified metabolites is over-represented in KEGG/SMPDB pathways, with a defensible background.

**Approach:** Map names/IDs to KEGG compounds and report mapping coverage; then run the hypergeometric test against an explicit assay-coverage background. Default to the **Local-Only ORA** function (no user data leaves the machine, background correction verified end to end); the MetaboAnalystR API path is the alternative.

```bash
# compounds.txt: one name/HMDB/KEGG ID per line (MSI level 1-2); writes kegg_ids.txt and prints mapping coverage
Rscript scripts/map_compounds.R compounds.txt hsa name kegg_ids.txt   # ID type: name | hmdb | kegg | pubchem
```

### Enrichment test

Run the hypergeometric test with the **Local-Only ORA** function (default) or the MetaboAnalystR API path; both are in `references/ora-local-and-api.md`.

## Reference Files

| File | Read when |
|---|---|
| `references/ora-local-and-api.md` | Running the ORA test on the mapped `kegg_ids`: Local-Only ORA (default, verified background correction) or the MetaboAnalystR API path (sends the list off-machine) |
| `references/mummichog-psea.md` | The input is a raw m/z peak table with no IDs: mummichog/PSEA code, ppm/ionization/p-cutoff setup |
| `references/fella-diffusion.md` | The goal is mechanism (enzymes/reactions linking the hits): FELLA graph build and `runDiffusion` |

## Per-Method Failure Modes

### Wrong background set (the silent controller)
- **Trigger:** ORA run with the full library ("all of KEGG"); mummichog run with only significant features as input.
- **Mechanism:** The background IS the null hypothesis made concrete. The KEGG-human library held ~3,373 compounds vs 286-1,110 actually measurable in real datasets; padding the denominator with undetectable compounds inflates every p-value. For mummichog, the permutation null samples from the input table, so a significant-only input pre-enriches the pool.
- **Symptom:** Many "significant" pathways; few survive once the background is the assay-specific metabolome (Wieder 2021: two of five datasets dropped to ZERO after FDR with the correct background).
- **Fix:** ORA -> Local-Only ORA (`references/ora-local-and-api.md`) with the measured-metabolome reference as `universe_kegg_ids` (or `SetMetabolomeFilter(mSet, TRUE)` after `Setup.KEGGReferenceMetabolome()` on the API path). Mummichog -> supply the entire feature table as `peaks.txt`. State the background in one sentence or the p-values are uninterpretable.

### Annotation laundering
- **Trigger:** ORA/MSEA run on MSI level-3 ("grey zone") tentative annotations as if they were level-1 confirmed.
- **Mechanism:** Enrichment cannot represent annotation confidence; a 4% misidentification rate already manufactures false pathways (Wieder 2021), and the error is not zero-mean because a wrong hub-adjacent compound flips a pathway alone.
- **Symptom:** Confident pathway claims downstream of unconfirmed IDs; results that do not replicate.
- **Fix:** Report the MSI level of the compounds driving the winning pathway; downgrade L3-driven claims to "consistent with." Consider metapone, which down-weights multiply-annotated features (weight inversely proportional to candidate count) instead of discarding the uncertainty.

### Hub-inflated topology / "impact"
- **Trigger:** Reporting MetaboAnalyst Pathway Impact as if it were an effect size.
- **Mechanism:** Impact = sum of relative-betweenness centrality of matched metabolites / sum over all pathway metabolites, computed inside an arbitrary isolated KEGG boundary without removing currency metabolites. A handful of cofactor-like hubs dominate betweenness (Tsouka & Masoodi 2023: L-alanine alone = ~95% of a pathway's total centrality). It is also sign-blind.
- **Symptom:** A pathway "lights up" with high impact because one promiscuous compound was hit by chance; the same hits give different impact in KEGG vs SMPDB.
- **Fix:** Treat impact as a visualization tiebreaker only. Read ORA and topology against each other; a pathway impact-driven by a single hub is a red flag, not a confirmation.

### Pool size is not flux
- **Trigger:** Reporting "pathway X is activated / upregulated" from concentration-based enrichment.
- **Mechanism:** A metabolomics measurement is a steady-state pool size (production minus consumption), not a rate. Pool and flux can move in opposite directions: sildenafil RAISES the cGMP pool while LOWERING flux through it (it inhibits the degrading phosphodiesterase). A falling substrate pool can mean the pathway is MORE active.
- **Symptom:** Causal/activity language ("upregulated pathway") drawn from a concentration snapshot.
- **Fix:** Downgrade to "members of pathway X co-varied with phenotype." Activity claims require stable-isotope-resolved metabolomics (SIRM / 13C metabolic flux analysis), which traces label incorporation over time; concentration-based enrichment generates a flux hypothesis, never tests one.

## Quantitative Thresholds

| Threshold | Value | Source / rationale |
|---|---|---|
| Mummichog query p-cutoff | ~0.2 (NOT 0.05) | The query must be large enough to score; vignette default is looser than 0.05. Document the value used (Li 2013; MetaboAnalystR vignette). |
| Empirical-compound RT window (v2) | ~`max(RT) * 0.02` seconds | Groups co-eluting features into one empirical compound; units are SECONDS (passing minutes mis-groups). |
| Mass tolerance (ppm) | instrument-specific (e.g. 5 ppm HRMS) | Loose ppm worsens multiple-m/z-matching inflation; set to the instrument's real accuracy. |
| FDR | < 0.05 (BH) | Standard, but secondary to a correct background -- with the right background, often zero pathways survive (Wieder 2021). |
| Pathway granularity caveat | -- | Pathway definition moves p by up to 9 orders of magnitude vs ~2 for multiple testing (Karp 2021); prefer cross-database consensus over one library. |
| Mapping coverage | report always | Enrichment computed over 12 of 400 features is a footnote, not a finding (Theme 3). |

## Common Errors

| Error / symptom | Cause | Solution |
|---|---|---|
| Everything is significant in mummichog | Input was significant features only, not R_all | Supply the entire feature table as `peaks.txt` |
| `could not find function "runPageRank"` | Wrong casing | FELLA function is `runPagerank` (lowercase r); method string is `'pagerank'` |
| PSEA maps to the wrong network silently | Wrong library string in `PerformPSEA` | Library encodes organism+network (`hsa_mfn`, `hsa_kegg`, ...); match the organism |
| Garbage candidate compounds | Wrong ionization mode | pos/neg use different adduct tables; set mode in `UpdateInstrumentParameters`; mixed data needs a per-feature mode column |
| Only TCA / amino-acid pathways enriched | Pathway dark matter | Xenobiotics, lipids, novel structures map to no pathway and are dropped; report coverage; consider ChemRICH (structure-based) |
| `'v2'` enrichment errors on RT | No RT column in input | `'v2'`/empirical compounds need RT; use `SetPeakFormat(mSet, 'mprt')` |
| `Error ... object 'current.msg' not found` on any validation failure (bad background, too few features, bad mapping) | `AddErrMsg()` crashes uncatchably outside the MetaboAnalyst web app | Pre-declare `current.msg <- character(0); err.vec <- character(0)` before any MetaboAnalystR call (see Version Compatibility); the real message then prints instead of crashing |

## References

- Li S, Park Y, Duraisingham S, Strobel FH, Khan N, Soltow QA, Jones DP, Pulendran B. 2013. Predicting network activity from high throughput metabolomics. *PLoS Comput Biol* 9(7):e1003123.
- Pang Z, Lu Y, Zhou G, Hui F, Xu L, Viau C, Spigelman AF, MacDonald PE, Wishart DS, Li S, Xia J. 2024. MetaboAnalyst 6.0: towards a unified platform for metabolomics data processing, analysis and interpretation. *Nucleic Acids Res* 52(W1):W398-W406.
- Xia J, Wishart DS. 2010. MSEA: a web-based tool to identify biologically meaningful patterns in quantitative metabolomic data. *Nucleic Acids Res* 38(W):W71-W77.
- Picart-Armada S, Fernandez-Albert F, Vinaixa M, Yanes O, Perera-Lluna A. 2018. FELLA: an R package to enrich metabolomics data. *BMC Bioinformatics* 19(1):538.
- Wieder C, Frainay C, Poupin N, Rodriguez-Mier P, Vinson F, Cooke J, Lai RPJ, Bundy JG, Jourdan F, Ebbels T. 2021. Pathway analysis in metabolomics: recommendations for the use of over-representation analysis. *PLOS Comput Biol* 17(9):e1009105.
- Wieder C, Bundy JG, Frainay C, Poupin N, Rodriguez-Mier P, Vinson F, Cooke J, Lai RPJ, Jourdan F, Ebbels TMD. 2022. Avoiding the misuse of pathway analysis tools in environmental metabolomics. *Environ Sci Technol* 56(20):14219-14222.
- Karp PD, Midford PE, Caspi R, Khodursky A. 2021. Pathway size matters: the influence of pathway granularity on over-representation (enrichment analysis) statistics. *BMC Genomics* 22:191.
- Tsouka S, Masoodi M. 2023. Metabolic pathway analysis: advantages and pitfalls for the functional interpretation of metabolomics and lipidomics data. *Biomolecules* 13(2):244.
- Tian L, Yu T. 2022. Metapone: a Bioconductor package for joint pathway testing for untargeted metabolomics data. *Bioinformatics* 38(14):3662-3669.
- Barupal DK, Fiehn O. 2017. Chemical Similarity Enrichment Analysis (ChemRICH) as alternative to biochemical pathway mapping for metabolomic datasets. *Sci Rep* 7:14567.
- Schymanski EL, Jeon J, Gulde R, Fenner K, Ruff M, Singer HP, Hollender J. 2014. Identifying small molecules via high resolution mass spectrometry: communicating confidence. *Environ Sci Technol* 48(4):2097-2098.

## Related Skills

- metabolomics/metabolite-annotation - Annotation confidence levels (MSI) that feed ORA/MSEA and set the interpretive ceiling
- metabolomics/statistical-analysis - Upstream differential testing that produces the significant compound or feature list
- metabolomics/isotope-tracing - Flux versus pool: enrichment infers activity from steady-state abundance, which isotope labeling can contradict
- pathway-analysis/go-enrichment - Gene-set over-representation concepts (the ORA analogue for genes)
- pathway-analysis/gsea - Ranked-list enrichment concepts (the MSEA analogue for genes)
- multi-omics-integration/mofa-integration - Joint gene+metabolite integration and its coverage-asymmetry traps
