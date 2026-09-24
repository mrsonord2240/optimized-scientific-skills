---
name: bio-proteomics-protein-inference
category: Data Analysis
description: Groups proteins from peptide identifications and controls protein-level FDR, framing inference as a chosen explanation (parsimony or a probability model) of underdetermined peptide evidence rather than a measurement. Reports protein GROUPS (proteins indistinguishable by observed peptides) with a leading protein, not flat lists. Covers shared-vs-unique peptides, indistinguishable/subsumable proteins, parsimony vs probabilistic (ProteinProphet, EPIFANY) vs razor inference, picked-protein and picked-group FDR, and why the two-peptide rule is wrong. Use when resolving which proteins are present from a peptide list, building protein groups, or estimating protein-level FDR. PSM/peptide FDR and search engines are peptide-identification; razor-vs-unique quant consequences are quantification; isoform/proteoform resolution is top-down and out of scope.
tool_type: mixed
primary_tool: pyOpenMS
license: MIT
author: GPTomics
---

## Version Compatibility

Reference examples tested with: pyOpenMS 3.5.0, pandas 2.2+, Percolator 3.09.0, Philosopher 5.1.0
Install: `pip install pyopenms pandas`; Percolator, OpenMS `Epifany` and Philosopher are separate CLI installs (Philosopher 5.1.0 ships the TPP-derived subcommands, so a separate TPP install is not needed).

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters
- CLI: `<tool> --help` (or `<tool> help <subcommand>`) and confirm the flag is in the listing BEFORE
  scripting it -- protein-inference flags get renamed and removed between releases (Percolator's Fido
  `--protein`/`-A` is the example; see the taxonomy)

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

In pyOpenMS the EPIFANY class is `BayesianProteinInferenceAlgorithm` (there is no `EpifanyAlgorithm`; the TOPP tool is `Epifany`). From 3.5, peptide IDs must be a `PeptideIdentificationList`, and `ProteinGroup.accessions` come back as bytes.

# Protein Inference -- A Chosen Explanation of Peptide Evidence, Reported as Groups

**"Tell me which proteins are present from my identified peptides"** -> Assign the observed peptides to a minimal or probability-weighted set of proteins, reported as groups of indistinguishable proteins with a leading accession -- because bottom-up MS measures peptides, and the protein set behind them is inferred, not observed.
- Python: `pyopenms.BasicProteinInferenceAlgorithm().run(peptide_ids, protein_ids)` for score-aggregation inference; set `greedy_group_resolution` to resolve shared peptides parsimony-style
- Python: `pyopenms.BayesianProteinInferenceAlgorithm` (TOPP tool `Epifany`) for Bayesian belief-propagation inference
- CLI: `percolator -f <target_decoy.fasta>` (`--picked-protein`) for picked-protein FDR straight off the `.pin` you already ran; `philosopher proteinprophet` + `philosopher filter --picked --razor` for the FragPipe/TPP route (EM apportionment) -- both judged by their OUTPUT FILES, never by their exit status

Scope: this skill OWNS peptide-to-protein grouping, the indistinguishable/subsumable distinction, the leading-protein convention, inference-method choice, and protein/protein-group FDR. PSM-level and peptide-level FDR plus the search engines that produce the peptide list -> peptide-identification. The quantitative fallout of razor vs unique peptides on protein abundance -> quantification. OUT OF SCOPE: resolving splice isoforms, single-AA variants, or PTM-defined proteoforms (bottom-up groups cannot separate them; that is top-down / proteoform work). Research use only: an inferred group is a chosen explanation, not a clinical finding -- route any patient-level "protein X is present" claim to a validated targeted assay (PRM/MRM on protein-unique peptides), never to a discovery protein list.

## The Single Most Important Modern Insight -- Protein Inference Is Underdetermined, So the Honest Unit Is a Group, Not a List

1. **The protein set is not uniquely recoverable from peptides, so a protein group -- not a flat protein list -- is the only honest reporting unit.** Many peptides are shared across paralogs, gene families, and isoforms, so distinct protein sets can explain the same peptide evidence equally well. The inference picks ONE explanation under an assumption (parsimony, or a probability model); proteins that the observed peptides cannot tell apart (indistinguishable) MUST be reported as one group with a designated leading protein. A flat list double-counts indistinguishable proteins and breaks target/decoy symmetry at the protein level, silently corrupting FDR.

2. **Protein FDR is its own estimation problem; estimate it at the protein level with PICKED FDR.** Controlling PSM-FDR at 1% does not give 1% protein-FDR: a deep run has many false PSMs in absolute terms, each can nucleate a one-hit-wonder false protein, and with NO protein-level estimate the protein list can be 10-30% false. That range is dataset- and threshold-dependent, not a constant: on this skill's own synthetic 113k-PSM benchmark, already filtered to 1% PSM FDR, the unfiltered protein list was 7.0% false (and the built-in estimator put it at 14.8%). Measure it on your own data instead of quoting the range. The classic (non-picked) protein-level target-decoy count errs the other way on large data: decoy proteins accumulate random hits faster than false targets, so it OVER-estimates FDR and discards real proteins. Savitski 2015 picked-protein FDR pairs each target protein with its decoy and keeps only the higher-scoring of the pair before counting, which removes that bias. Picked-group FDR (The, Samaras, Kuster & Wilhelm 2022) applies picking to protein groups, and also shows that naive Occam/parsimony-style grouping can be anticonservative on large data. Subsumable proteins must be removed before groups are counted, or they inflate the target list (synthetic test: 8.6% true FDP at nominal 1% without resolution, 0.9% with it).

3. **Do not use the two-peptide rule -- it discards real proteins, and it is not an FDR control.** Requiring >=2 peptides per protein (Gupta & Pevzner 2009, "A strike against the two-peptide rule") throws away legitimate low-abundance single-peptide IDs. Its effect on protein FDR depends on the PSM threshold: Gupta & Pevzner found it raised FDR, while under a strict 1% PSM pre-filter it can lower it at the cost of many true proteins. Replace the blanket rule with: control protein-level (picked) FDR, then judge single-peptide IDs by their score, not their peptide count.

## Vocabulary the Rest of This Depends On

- Shared (degenerate) peptide: maps to >1 protein in the searched database. Cannot, alone, distinguish which protein is present.
- Unique peptide: maps to exactly one protein -- the only direct evidence for a specific protein. "Unique" is DATABASE-RELATIVE: a peptide unique against SwissProt may be shared against TrEMBL+isoforms+contaminants. Always document the exact database (isoforms, contaminants, decoys included).
- Indistinguishable proteins: explained by the SAME set of observed peptides -> one group, never two confident IDs.
- Subset / subsumable protein: its observed peptides are a subset of another protein's -> parsimony drops it (the larger protein explains everything it would).
- Leading / representative protein: the group's reported accession. Convention: most peptides, then highest score, then SwissProt canonical over TrEMBL. Downstream tables key on this accession but must retain group membership -- "protein P12345" usually means "the group led by P12345".
- Protein group vs proteoform: a group is an inference artifact (proteins lumped because peptides cannot separate them); a proteoform is a real molecular species (one gene product with a specific sequence + PTM + cleavage state). Bottom-up groups DO NOT resolve proteoforms -- claiming "isoform X present" from a shared-peptide group is overreach.

## Tool Taxonomy

| Tool / method | Citation | Mechanism / role | When |
|---------------|----------|------------------|------|
| Parsimony (Occam) | -- | Greedy minimal protein set explaining all peptides | Fast default; ties broken arbitrarily; can be anticonservative for group FDR on large data (The 2022) |
| Score aggregation (OpenMS `BasicProteinInferenceAlgorithm`) | -- | Aggregates peptide scores for EVERY protein; `greedy_group_resolution` adds razor-style resolution | Not parsimony by default: subsumable and shared-only proteins stay as groups unless `greedy_group_resolution='true'` |
| ProteinProphet | Nesvizhskii 2003 | EM APPORTIONS shared peptides across candidate proteins, weighted by other evidence | TPP / FragPipe pipelines; the classic probabilistic standard |
| EPIFANY | Pfeuffer 2020 | Bayesian network over the peptide-protein graph, loopy belief propagation + convolution trees | When you want calibrated posteriors rather than a parsimony set. It is NOT automatically better calibrated at the group level: on the reference idXML it passed 583 groups at 1% picked FDR with 3.43% true FDP, against 556 and 0.36% for Basic + `greedy_group_resolution` |
| Fido | Serang 2010 | Bayesian generative model; was Percolator's `--protein` / `-A` option | GONE from the Percolator CLI -- absent from 3.09.0 `--help`; `--protein` exits 1 (see Common Errors). Only "Not available for Fido" remarks on `--protein-report-fragments`/`-duplicates` survive. Do not script it |
| Percolator picked-protein (`-f` / `--picked-protein <fasta>`) | Savitski 2015 | In-silico digest of the FASTA -> protein grouping -> fragment/duplicate elimination -> picking, on the `.pin` Percolator already read | The way Percolator 3.09 does protein inference, and it IS the picked method recommended here -- no extra tool needed after Percolator. Default output is one representative per row, NOT a group list (partners are eliminated; see `references/percolator-philosopher-cli.md`) |
| Razor peptide | -- | Shared peptide assigned winner-take-all to the group with most evidence (MaxQuant) | MaxQuant default; ID-fine but distorts QUANT (route to quantification) |
| Picked-protein FDR | Savitski 2015 | Pair target with its decoy, keep the higher-scoring of the pair, then count decoys | Protein-level FDR on any non-trivial dataset |
| Picked-group FDR | The 2022 | Picking applied at the protein-GROUP level | When the inference unit is the group (the correct unit on deep data) |
| All-proteins / inclusive | -- | Report every protein any peptide could come from | Almost never; massive false-positive protein inflation |

## Decision Tree by Scenario

| Scenario | Recommended | Why |
|----------|-------------|-----|
| Standard DDA run, OpenMS-based pipeline | `BasicProteinInferenceAlgorithm` with `greedy_group_resolution='true'` + picked-group FDR; `BayesianProteinInferenceAlgorithm` (EPIFANY) when you want calibrated posteriors | Both are group-FDR aware, but resolution is what controls FDP: on the same idXML, Basic+greedy gave 0.36% true FDP and 0 subsumable groups passing, EPIFANY called with `greedy_group_resolution=False` gave 3.43% and 3 (code: `references/pyopenms-basic-inference.md`, `references/epifany-bayesian-inference.md`) |
| Already running Percolator after the search | `percolator -f target_decoy.fasta -P DECOY_ -l prot.target.tsv -L prot.decoy.tsv ... search.pin` | Percolator's own picked-protein route; adds no tool and computes the method recommended below. Add `--protein-report-duplicates --protein-report-fragments` if the report must list group members (see `references/percolator-philosopher-cli.md`) |
| MaxQuant output (`proteinGroups.txt`) | Parse groups as-is; drop `Reverse` (`+`, accessions `REV__`), `Potential contaminant` and `Only identified by site` rows; quantify on UNIQUE peptides | Groups already inferred; `Protein IDs` = all members, `Majority protein IDs` = members with at least half the group's peptides, first = leading; `Unique peptides` = unique to the GROUP, not to one protein |
| FragPipe / TPP pipeline | `philosopher peptideprophet` -> `proteinprophet` -> `filter --picked --razor` -> `report` (commands in `references/percolator-philosopher-cli.md`), then COUNT THE ROWS in `protein.tsv` | Native EM apportionment + 2-level FDR; `filter` exits 0 and prints `Converged to 0.00 % FDR` when it read nothing, so the row count is the only real result |
| Deep dataset (many thousands of proteins) | Picked-GROUP FDR on resolved groups (`references/picked-group-fdr.md`) | No protein-level FDR leaves a 7-30% false protein list (measure it; the figure depends on depth and the PSM threshold); the non-picked count over-estimates FDR and loses proteins |
| Sensitive differential abundance downstream | Quantify on unique peptides only -> quantification | Razor assignment can flip between conditions and fake DE |
| Want isoform-level answers | Stop -- route to top-down / proteoform methods | Bottom-up groups cannot resolve proteoforms |
| Few PSMs (single-protein pulldown) | Report evidence, do not trust a "0% protein FDR" | Target-decoy FDR is meaningless at tiny counts |

Default when uncertain: run `BasicProteinInferenceAlgorithm` with `annotate_indistinguishable_groups` and `greedy_group_resolution` on, report protein GROUPS with a leading accession, and control protein-GROUP FDR with picked-group FDR at 1% (`references/pyopenms-basic-inference.md`). Do NOT impose a two-peptide rule.

**Input contract:** protein FDR needs decoys, so the upstream PSM/peptide filter must KEEP the decoy hits that pass it. Know the score orientation: PEP is lower-better, posterior probability and the group probability written by these algorithms are higher-better.

**Output contract (CLI):** judge every command-line inference/FDR step by the file it was supposed to write, NEVER by its exit status. A protein-FDR tool that read zero PSMs reports the emptiest possible result as a triumph (`philosopher filter` and `philosopher peptideprophet` both exit 0 on it; reproduced on Philosopher 5.1.0 with real Comet pepXML -- see Common Errors). After each step assert non-empty: `protein.tsv` has more than its header line, the `prot.xml` exists, the peptide count Percolator prints is greater than zero. A "0.00 % FDR" or a "0% protein FDR" is a failed run until a row count says otherwise.

## Reference Files

Each file holds one method's goal, approach and runnable code. SKILL.md is enough to choose a method; read the file the decision tree points at before running it.

| File | Read when |
|------|-----------|
| `references/pyopenms-basic-inference.md` | Default OpenMS route: `BasicProteinInferenceAlgorithm` + `greedy_group_resolution` + picked-group FDR, reported as group records |
| `references/epifany-bayesian-inference.md` | Calibrated posteriors wanted from PSMs that already carry PEPs (`BayesianProteinInferenceAlgorithm`) |
| `references/picked-group-fdr.md` | Estimating protein-group FDR by picking; decoy-prefix errors |
| `references/percolator-philosopher-cli.md` | Percolator `-f` picked-protein route or the Philosopher/FragPipe route, with the output check after every step |

## Per-Method Failure Modes

### Naive (non-picked) protein/group FDR
**Trigger:** (a) No protein-level FDR at all ("1% PSM FDR is enough"); (b) the classic protein-level `decoys/targets` count without picking on a deep dataset; (c) picked FDR on unresolved groups.
**Mechanism:** (a) false PSMs nucleate one-hit-wonder false proteins that no protein-level estimate catches; (b) decoy proteins keep accumulating random hits while true targets saturate, so decoys are over-counted; (c) subsumable and shared-only proteins count as extra target groups.
**Symptom:** (a) protein list 7-30% false depending on depth and PSM threshold (7.0% observed on the 113k-PSM synthetic benchmark at 1% PSM FDR); (b) FDR over-estimated, real proteins lost (conservative, not anticonservative); (c) nominal 1% with several-fold higher true FDP.
**Fix:** Picked-protein FDR (Savitski 2015) or picked-group FDR (The 2022) on resolved groups; validate with an entrapment search, which needs a proteome ABSENT from the sample (append e.g. Arabidopsis or a shuffled second proteome to the search database) -- a multi-species benchmark such as a HYE human/yeast/E. coli mix does not provide one, because every species in it is truly present.

### Two-peptide rule
**Trigger:** Filtering to proteins with >=2 (unique) peptides "for confidence".
**Mechanism:** The rule deletes real low-abundance single-peptide proteins; its FDR effect depends on the PSM threshold (Gupta & Pevzner found it raised FDR; after a strict 1% PSM filter it can lower it while still deleting hundreds of true proteins).
**Symptom:** Fewer proteins than picked FDR at the same nominal cutoff, with no calibrated error rate.
**Fix:** Drop the rule; control picked protein-level FDR and score single-peptide IDs individually.

### Razor-peptide quantification
**Trigger:** Quantifying on MaxQuant's default unique+razor peptides for a sensitive comparison.
**Mechanism:** A shared peptide's full intensity is credited to one group; that razor assignment can flip between conditions when peptide counts shift, so a protein's quantity changes for inference reasons, not biology.
**Symptom:** Spurious differential abundance concentrated on proteins sharing peptides with paralogs.
**Fix:** Quantify on unique peptides only for sensitive comparisons -> quantification.

### Parsimony tie-breaking
**Trigger:** Multiple minimal protein sets explain the peptides equally well.
**Mechanism:** Greedy parsimony breaks ties arbitrarily; minimality is a heuristic, not truth, and a real protein with only shared peptides is silently dropped.
**Symptom:** Reported lead protein differs run-to-run or pipeline-to-pipeline on the same data.
**Fix:** Prefer a probabilistic method (EPIFANY/ProteinProphet) that apportions shared evidence; retain group membership.

### Proteoform overreach
**Trigger:** Reporting "isoform X is present" from a group whose evidence is shared peptides.
**Mechanism:** Splice isoforms, variants, and PTM forms collapse into groups in bottom-up data; the group cannot separate them.
**Symptom:** Isoform-specific claim with no isoform-unique peptide behind it.
**Fix:** Require an isoform-unique peptide for any isoform claim, or use top-down / proteoform methods.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| Protein / protein-group FDR 1% (sometimes 5% for discovery) | community standard | SEPARATE estimation from PSM FDR; never assume 1% PSM implies 1% protein |
| Picked FDR (target/decoy pairing) | Savitski 2015; The 2022 | Removes target/decoy asymmetry; dataset-size-independent, unlike naive decoy/target |
| Decoy:target ratio 1:1 | community standard | Standard null; unequal ratios require formula correction |
| Min PSMs for trustworthy protein FDR | hundreds+ | Below ~100s of items decoy counts are too noisy; "0% FDR" from zero decoys is luck, not control |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| `AttributeError: module 'pyopenms' has no attribute 'EpifanyAlgorithm'` | The pyOpenMS class is named differently; the R `ProteinInference::infer_proteins` could not be confirmed to exist | Use `pyopenms.BayesianProteinInferenceAlgorithm`; use pyOpenMS, not an unverified R package |
| `Exception: can not handle type of (..., [], [])` on `IdXMLFile().load` | pyOpenMS 3.5+ needs a `PeptideIdentificationList` | `peptide_ids = PeptideIdentificationList()` |
| `TypeError: a bytes-like object is required, not 'str'` on group accessions | `ProteinGroup.accessions` are bytes | `[a.decode() for a in group.accessions]` |
| `IndexError: invalid unordered_map<K, T> key` from `applyPickedProteinFDR` (or `ValueError: no decoy groups with prefix ...` from `scripts/picked_group_fdr.py`) | The decoy prefix does not match the one in the file | Pass the tool's prefix: `DECOY_` OpenMS, `rev_` Philosopher/FragPipe, `REV__` MaxQuant. The error means "wrong prefix", not "corrupt input" |
| `ERROR: the option --protein is invalid.` (exit 1) from Percolator | Fido's `--protein`/`-A` was removed; it is not in 3.09.0 `--help` | `percolator -f target_decoy.fasta -P <decoy prefix> -l prot.target.tsv ...` (`--picked-protein`) |
| `philosopher filter` exits 0, logs `Converged to 0.00 % FDR with 0 PSMs` and `proteins=0`, writes header-only or no tables | It read zero PSMs from the pepXML; the exit status reports the empty result as success | Never trust the exit code -- assert `protein.tsv` has more than a header line. Then fix the input: check `Database search results ions=/peptides=/psms=` in the log, that `--tag` matches the FASTA prefix (default is `rev_`), and that `philosopher database --annotate` ran first |
| `philosopher peptideprophet` exits 0 but `proteinprophet` then says `did not find any PeptideProphet results` | PeptideProphet modelled nothing and still wrote an `interact-*.pep.xml` full of `spectrum_query` elements with no probabilities | `grep -c peptideprophet_result interact-*.pep.xml` immediately after the step; 0 means the run failed |
| Indistinguishable proteins reported as separate IDs | Flat protein list instead of groups | Enable `annotate_indistinguishable_groups`; report groups with a leading protein |

## References

- Nesvizhskii, A.I., Keller, A., Kolker, E. & Aebersold, R. (2003). A statistical model for identifying proteins by tandem mass spectrometry. *Analytical Chemistry* 75(17):4646-4658.
- Gupta, N. & Pevzner, P.A. (2009). False discovery rates of protein identifications: a strike against the two-peptide rule. *Journal of Proteome Research* 8(9):4173-4181.
- Serang, O., MacCoss, M.J. & Noble, W.S. (2010). Efficient marginalization to compute protein posterior probabilities from shotgun mass spectrometry data. *Journal of Proteome Research* 9(10):5346-5357. (Fido; its Percolator front-end no longer exists as of 3.09.0.)
- Savitski, M.M., Wilhelm, M., Hahne, H., Kuster, B. & Bantscheff, M. (2015). A scalable approach for protein false discovery rate estimation in large proteomic data sets. *Molecular & Cellular Proteomics* 14(9):2394-2404.
- The, M., Tasnim, A. & Kall, L. (2016). How to talk about protein-level false discovery rates in shotgun proteomics. *Proteomics* 16(18):2461-2469.
- The, M., Samaras, P., Kuster, B. & Wilhelm, M. (2022). Reanalysis of ProteomicsDB using an accurate, sensitive, and scalable false discovery rate estimation approach for protein groups. *Molecular & Cellular Proteomics* 21(12):100437.
- Pfeuffer, J., Sachsenberg, T., Dijkstra, T.M.H., Serang, O., Reinert, K. & Kohlbacher, O. (2020). EPIFANY: a method for efficient high-confidence protein inference. *Journal of Proteome Research* 19(3):1060-1072.

## Related Skills

- peptide-identification - Produces the FDR-filtered peptide list that feeds inference and shares the target-decoy machinery
- quantification - Consumes inferred groups; razor-vs-unique peptide choice lives here
- data-import - Loads idXML/mzML identification files
- database-access/uniprot-access - Canonical-vs-isoform databases drive uniqueness and the leading-protein convention
