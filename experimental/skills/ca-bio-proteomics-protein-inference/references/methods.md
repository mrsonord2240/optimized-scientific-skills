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
