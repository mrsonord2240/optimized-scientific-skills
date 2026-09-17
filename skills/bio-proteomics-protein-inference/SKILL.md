---
name: bio-proteomics-protein-inference
description: Groups proteins from peptide identifications and controls protein-level FDR, framing inference as a chosen explanation (parsimony or a probability model) of underdetermined peptide evidence rather than a measurement. Reports protein GROUPS (proteins indistinguishable by observed peptides) with a leading protein, not flat lists. Covers shared-vs-unique peptides, indistinguishable/subsumable proteins, parsimony vs probabilistic (ProteinProphet, EPIFANY) vs razor inference, picked-protein and picked-group FDR, and why the two-peptide rule is wrong. Use when resolving which proteins are present from a peptide list, building protein groups, or estimating protein-level FDR. PSM/peptide FDR and search engines are peptide-identification; razor-vs-unique quant consequences are quantification; isoform/proteoform resolution is top-down and out of scope.
tool_type: mixed
primary_tool: pyOpenMS
license: MIT
---

## Version Compatibility

Reference examples tested with: pyOpenMS 3.5.0, pandas 2.2+, Percolator 3.09.0, Philosopher 5.1.0

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters
- CLI: `<tool> --help` (or `<tool> help <subcommand>`) and confirm the flag is in the listing BEFORE
  scripting it -- protein-inference flags get renamed and removed between releases (Percolator's Fido
  option `--protein`/`-A` is gone in 3.09.0 and exits 1 with `ERROR: the option --protein is invalid.`)

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
| Fido | Serang 2010 | Bayesian generative model; was Percolator's `--protein` / `-A` option | GONE from the Percolator CLI -- absent from 3.09.0 `--help`, and `--protein` exits 1 with `ERROR: the option --protein is invalid.` Only "Not available for Fido" remarks on `--protein-report-fragments`/`-duplicates` survive. Do not script it |
| Percolator picked-protein (`-f` / `--picked-protein <fasta>`) | Savitski 2015 | In-silico digest of the FASTA -> protein grouping -> fragment/duplicate elimination -> picking, on the `.pin` Percolator already read | The way Percolator 3.09 does protein inference, and it IS the picked method recommended here -- no extra tool needed after Percolator |
| Razor peptide | -- | Shared peptide assigned winner-take-all to the group with most evidence (MaxQuant) | MaxQuant default; ID-fine but distorts QUANT (route to quantification) |
| Picked-protein FDR | Savitski 2015 | Pair target with its decoy, keep the higher-scoring of the pair, then count decoys | Protein-level FDR on any non-trivial dataset |
| Picked-group FDR | The 2022 | Picking applied at the protein-GROUP level | When the inference unit is the group (the correct unit on deep data) |
| All-proteins / inclusive | -- | Report every protein any peptide could come from | Almost never; massive false-positive protein inflation |

## Decision Tree by Scenario

| Scenario | Recommended | Why |
|----------|-------------|-----|
| Standard DDA run, OpenMS-based pipeline | `BasicProteinInferenceAlgorithm` with `greedy_group_resolution='true'` + picked-group FDR; `BayesianProteinInferenceAlgorithm` (EPIFANY) when you want calibrated posteriors | Both are group-FDR aware, but resolution is what controls FDP: on the same idXML, Basic+greedy gave 0.36% true FDP and 0 subsumable groups passing, EPIFANY called with `greedy_group_resolution=False` gave 3.43% and 3 |
| Already running Percolator after the search | `percolator -f target_decoy.fasta -P DECOY_ -l prot.target.tsv -L prot.decoy.tsv ... search.pin` | Percolator's own picked-protein route; adds no tool and computes the method recommended below |
| MaxQuant output (`proteinGroups.txt`) | Parse groups as-is; drop `Reverse` (`+`, accessions `REV__`), `Potential contaminant` and `Only identified by site` rows; quantify on UNIQUE peptides | Groups already inferred; `Protein IDs` = all members, `Majority protein IDs` = members with at least half the group's peptides, first = leading; `Unique peptides` = unique to the GROUP, not to one protein |
| FragPipe / TPP pipeline | `philosopher peptideprophet` -> `proteinprophet` -> `filter --picked --razor` -> `report` (commands below), then COUNT THE ROWS in `protein.tsv` | Native EM apportionment + 2-level FDR; `filter` exits 0 and prints `Converged to 0.00 % FDR` when it read nothing, so the row count is the only real result |
| Deep dataset (many thousands of proteins) | Picked-GROUP FDR on resolved groups | No protein-level FDR leaves a 7-30% false protein list (measure it; the figure depends on depth and the PSM threshold); the non-picked count over-estimates FDR and loses proteins |
| Sensitive differential abundance downstream | Quantify on unique peptides only -> quantification | Razor assignment can flip between conditions and fake DE |
| Want isoform-level answers | Stop -- route to top-down / proteoform methods | Bottom-up groups cannot resolve proteoforms |
| Few PSMs (single-protein pulldown) | Report evidence, do not trust a "0% protein FDR" | Target-decoy FDR is meaningless at tiny counts |

Default when uncertain: run `BasicProteinInferenceAlgorithm` with `annotate_indistinguishable_groups` and `greedy_group_resolution` on, report protein GROUPS with a leading accession, and control protein-GROUP FDR with picked-group FDR at 1%. Do NOT impose a two-peptide rule.

**Input contract:** protein FDR needs decoys, so the upstream PSM/peptide filter must KEEP the decoy hits that pass it. Know the score orientation: PEP is lower-better, posterior probability and the group probability written by these algorithms are higher-better.

**Output contract (CLI):** judge every command-line inference/FDR step by the file it was supposed to write, NEVER by its exit status. A protein-FDR tool that read zero PSMs reports the emptiest possible result as a triumph -- `philosopher filter` exits 0 and prints `Converged to 0.00 % FDR with 0 PSMs` and `Final report ... proteins=0`, and `philosopher peptideprophet` exits 0 after writing an `interact-*.pep.xml` with every `spectrum_query` copied through and zero `peptideprophet_result` elements (both reproduced on Philosopher 5.1.0 with real Comet pepXML). After each step assert non-empty: `protein.tsv` has more than its header line, the `prot.xml` exists, the peptide count Percolator prints is greater than zero. A "0.00 % FDR" or a "0% protein FDR" is a failed run until a row count says otherwise.

### Group Proteins with pyOpenMS (aggregation + greedy resolution)

**Goal:** Turn an FDR-filtered peptide identification list into protein groups with a leading protein, resolving shared-peptide ambiguity.

**Approach:** Load the idXML from peptide identification, run `BasicProteinInferenceAlgorithm` (score aggregation per protein) with indistinguishable-group annotation AND greedy group resolution on -- without resolution, subsumable and shared-only proteins stay as their own groups -- then read the groups off the protein identification run and apply picked protein-group FDR with the built-in.

```python
from pyopenms import (IdXMLFile, BasicProteinInferenceAlgorithm, PeptideIdentificationList,
                      FalseDiscoveryRate, String)

protein_ids = []
peptide_ids = PeptideIdentificationList()   # pyOpenMS 3.5+: a plain [] fails
# protein_ids is FIRST in both load() and store() for IdXMLFile
IdXMLFile().load('peptides_1pct_fdr.idXML', protein_ids, peptide_ids)

inference = BasicProteinInferenceAlgorithm()
params = inference.getParameters()
# annotate_indistinguishable_groups reports indistinguishable proteins as ONE group
params.setValue('annotate_indistinguishable_groups', 'true')
# greedy_group_resolution assigns shared peptides to the best group, so subsumable proteins drop out
params.setValue('greedy_group_resolution', 'true')
inference.setParameters(params)
inference.run(peptide_ids, protein_ids)

# picked protein-group FDR: (decoy string, is prefix, groups too); group.probability becomes a q-value
FalseDiscoveryRate().applyPickedProteinFDR(protein_ids[0], String('DECOY_'), True, True)

for group in protein_ids[0].getIndistinguishableProteins():
    accs = [a.decode() for a in group.accessions]   # bytes; pyOpenMS sorts them alphabetically
    if all(a.startswith('DECOY_') for a in accs) or group.probability > 0.01:
        continue
    # members share the same evidence; choose the lead explicitly: canonical (no -N isoform suffix) first
    leading = sorted(accs, key=lambda a: ('-' in a, a))[0]
    print(leading, group.probability, accs)
```

### Bayesian Inference + Group FDR with EPIFANY

**Goal:** Assign calibrated protein/group posteriors and control protein-group FDR with a probability model rather than greedy parsimony.

**Approach:** EPIFANY consumes idXML whose PSMs already carry posterior error probabilities (from Percolator or IDPosteriorErrorProbability), then propagates belief over the peptide-protein graph. The TOPP tool is `Epifany`; the pyOpenMS class is `BayesianProteinInferenceAlgorithm`.

```python
from pyopenms import IdXMLFile, BayesianProteinInferenceAlgorithm, PeptideIdentificationList

protein_ids = []
peptide_ids = PeptideIdentificationList()
IdXMLFile().load('peptides_with_pep.idXML', protein_ids, peptide_ids)

algo = BayesianProteinInferenceAlgorithm()
# EPIFANY expects PSM posteriors as input. The third POSITIONAL argument is
# greedy_group_resolution. Unlike BasicProteinInferenceAlgorithm's parameter of the same
# name, flipping it here did not change the result on the reference idXML (583 groups at 1%
# picked FDR, 3.43% true FDP with either value) -- so do NOT assume it removed subsumable
# proteins. Check the surviving groups, or use Basic + greedy when FDP is what you care about.
algo.inferPosteriorProbabilities(protein_ids, peptide_ids, False)

for group in protein_ids[0].getIndistinguishableProteins():
    print([a.decode() for a in group.accessions], group.probability)   # higher = more likely present
```

### Picked Protein-Group FDR

**Goal:** Estimate protein-group FDR without the inflation that the reused PSM formula causes on large data.

**Approach:** For each target group, find its decoy counterpart (same accessions with the decoy prefix); keep only the higher-scoring member of each target/decoy PAIR; rank the picked set and count decoys as the FDR estimate. For idXML, prefer the built-in `FalseDiscoveryRate().applyPickedProteinFDR(prot_id, String(prefix), True, True)` shown above; for group-level work on large data, the kusterlab `picked_group_fdr` package implements The 2022. The sketch below pairs by the exact accession set, so decoy groups whose membership differs from their target's stay unpaired and are counted unpicked. The decoy prefix is tool-specific (`DECOY_` OpenMS, `rev_` FragPipe/Philosopher's `--tag` default, `REV__` MaxQuant) and must be passed explicitly. A wrong prefix does not pass silently in either implementation: `applyPickedProteinFDR` raises `IndexError: invalid unordered_map<K, T> key`, and the sketch below raises `ValueError: no decoy groups with prefix ...`. Treat both as "the prefix is wrong", not as a corrupt input file.

```python
def picked_group_fdr(groups, decoy_prefix, min_decoys=10):
    # groups: list of dicts with 'accessions' (str), 'score' (higher = better), 'is_decoy'
    n_decoy = sum(g['is_decoy'] for g in groups)
    if n_decoy == 0 or not any(a.startswith(decoy_prefix) for g in groups for a in g['accessions']):
        raise ValueError(f'no decoy groups with prefix {decoy_prefix!r}; keep decoys upstream or fix the prefix')
    if n_decoy < min_decoys:
        print(f'WARNING: only {n_decoy} decoy groups; the protein FDR estimate is not meaningful')
    by_base = {}
    for g in groups:
        base = frozenset(a.replace(decoy_prefix, '') for a in g['accessions'])
        # keep only the higher-scoring of the target/decoy pair (the 'pick')
        if base not in by_base or g['score'] > by_base[base]['score']:
            by_base[base] = g
    picked = sorted(by_base.values(), key=lambda g: g['score'], reverse=True)

    targets = decoys = 0
    for g in picked:
        if g['is_decoy']:
            decoys += 1
        else:
            targets += 1
        g['fdr'] = decoys / targets if targets else 1.0
    running_min = 1.0
    for g in reversed(picked):  # monotone q-values from the bottom up
        running_min = min(running_min, g['fdr'])
        g['qvalue'] = running_min
    return [g for g in picked if not g['is_decoy'] and g['qvalue'] <= 0.01]
```

### Protein Groups from the CLI: Percolator and Philosopher

**Goal:** Get groups and protein-level FDR out of a search you are already post-processing on the command line, and recognise a failed run instead of reporting its "0.00 % FDR".

**Approach:** Percolator does picked-protein inference itself from the `.pin` -- `-f/--picked-protein` takes the TARGET+DECOY FASTA, in-silico-digests it to build protein groups, eliminates fragment (subsumable) and duplicate proteins, then picks. The FragPipe/TPP route chains PeptideProphet -> ProteinProphet -> `philosopher filter`. Every step below is checked on the file it wrote (see the Output contract above).

```bash
# --- Percolator 3.09.0: picked-protein FDR, no extra tool after Percolator ---
# Fido (the old `--protein` / `-A`) is GONE in 3.09.0; `-f` is the only protein route.
# -P is the decoy prefix in the FASTA; -z must match the search enzyme (default trypsin).
percolator -f target_decoy.fasta -P DECOY_ -z trypsin \
           -l prot.target.tsv -L prot.decoy.tsv \
           -r pep.target.tsv  -B pep.decoy.tsv \
           -m psm.target.tsv  -M psm.decoy.tsv -S 1 search.pin
# OUTPUT CHECK. prot.target.tsv columns: ProteinId, ProteinGroupId, q-value,
# posterior_error_prob, peptideIds -- one row per group representative, already picked.
awk -F'\t' 'NR>1 && $3<=0.01' prot.target.tsv | wc -l    # groups at 1% protein-group FDR

# --- Philosopher 5.1.0: the FragPipe / TPP route ---
philosopher workspace --init                              # once per output folder
philosopher database --annotate target_decoy.fasta --prefix DECOY_
philosopher peptideprophet --database target_decoy.fasta --ppm --accmass \
                           --nonparam --decoy DECOY_ --decoyprobs search.pep.xml
# OUTPUT CHECK: peptideprophet exits 0 even when it modelled nothing, copying every
# spectrum_query through with no probability attached. Count the results, not the exit code.
grep -c peptideprophet_result interact-search.pep.xml || { echo 'PeptideProphet modelled 0 PSMs'; exit 1; }
philosopher proteinprophet --maxppmdiff 2000000 interact-search.pep.xml
test -s interact.prot.xml || { echo 'ProteinProphet wrote no prot.xml'; exit 1; }
# --tag defaults to 'rev_', NOT the OpenMS/Comet 'DECOY_'; --razor is silently ignored
# unless --protxml supplies inference data.
philosopher filter --psm 0.01 --pep 0.01 --prot 0.01 --tag DECOY_ --picked --razor \
                   --pepxml interact-search.pep.xml --protxml interact.prot.xml
philosopher report
# OUTPUT CHECK: `filter` exits 0 and prints "Converged to 0.00 % FDR with 0 PSMs" on an
# empty read, and writes header-only tables or none at all.
[ -s protein.tsv ] && [ "$(wc -l < protein.tsv)" -gt 1 ] || { echo 'filter converged on an EMPTY result'; exit 1; }
```

## Per-Method Failure Modes

### Naive (non-picked) protein/group FDR
**Trigger:** (a) No protein-level FDR at all ("1% PSM FDR is enough"); (b) the classic protein-level `decoys/targets` count without picking on a deep dataset; (c) picked FDR on unresolved groups.
**Mechanism:** (a) false PSMs nucleate one-hit-wonder false proteins that no protein-level estimate catches; (b) decoy proteins keep accumulating random hits while true targets saturate, so decoys are over-counted; (c) subsumable and shared-only proteins count as extra target groups.
**Symptom:** (a) protein list 7-30% false depending on depth and PSM threshold (7.0% observed on the 113k-PSM synthetic benchmark at 1% PSM FDR); (b) FDR over-estimated, real proteins lost (conservative, not anticonservative); (c) nominal 1% with several-fold higher true FDP.
**Fix:** Picked-protein FDR (Savitski 2015) or picked-group FDR (The 2022) on resolved groups; validate with a two-species or entrapment search.

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
| Two-peptide rule | DO NOT USE (Gupta & Pevzner 2009) | Drops real proteins; FDR effect depends on the PSM threshold; replaced by picked FDR + per-ID score |
| Single-peptide IDs | judge by score, not count | A high-confidence unique peptide can be a legitimate ID |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Protein FDR much higher than nominal on deep data | No protein-level estimate, or groups not resolved (Basic without `greedy_group_resolution`) | Resolve groups, then picked-protein or picked-group FDR |
| Real low-abundance proteins missing | Two-peptide rule applied, or non-picked protein FDR (conservative) | Remove the rule; control picked FDR |
| `AttributeError: module 'pyopenms' has no attribute 'EpifanyAlgorithm'` | The pyOpenMS class is named differently; the R `ProteinInference::infer_proteins` could not be confirmed to exist | Use `pyopenms.BayesianProteinInferenceAlgorithm`; use pyOpenMS, not an unverified R package |
| `Exception: can not handle type of (..., [], [])` on `IdXMLFile().load` | pyOpenMS 3.5+ needs a `PeptideIdentificationList` | `peptide_ids = PeptideIdentificationList()` |
| `TypeError: a bytes-like object is required, not 'str'` on group accessions | `ProteinGroup.accessions` are bytes | `[a.decode() for a in group.accessions]` |
| `IndexError: invalid unordered_map<K, T> key` from `applyPickedProteinFDR` (or `ValueError: no decoy groups with prefix ...` from the sketch) | The decoy prefix does not match the one in the file | Pass the tool's prefix: `DECOY_` OpenMS, `rev_` Philosopher/FragPipe, `REV__` MaxQuant. The error means "wrong prefix", not "corrupt input" |
| `ERROR: the option --protein is invalid.` (exit 1) from Percolator | Fido's `--protein`/`-A` was removed; it is not in 3.09.0 `--help` | `percolator -f target_decoy.fasta -P <decoy prefix> -l prot.target.tsv ...` (`--picked-protein`) |
| `philosopher filter` exits 0, logs `Converged to 0.00 % FDR with 0 PSMs` and `proteins=0`, writes header-only or no tables | It read zero PSMs from the pepXML; the exit status reports the empty result as success | Never trust the exit code -- assert `protein.tsv` has more than a header line. Then fix the input: check `Database search results ions=/peptides=/psms=` in the log, that `--tag` matches the FASTA prefix (default is `rev_`), and that `philosopher database --annotate` ran first |
| `philosopher peptideprophet` exits 0 but `proteinprophet` then says `did not find any PeptideProphet results` | PeptideProphet modelled nothing and still wrote an `interact-*.pep.xml` full of `spectrum_query` elements with no probabilities | `grep -c peptideprophet_result interact-*.pep.xml` immediately after the step; 0 means the run failed |
| Indistinguishable proteins reported as separate IDs | Flat protein list instead of groups | Enable `annotate_indistinguishable_groups`; report groups with a leading protein |
| Spurious DE on paralog-sharing proteins | Razor-peptide quant flipped between conditions | Quantify on unique peptides -> quantification |
| "Unique" peptide count changed when DB changed | Uniqueness is database-relative | Fix and document the database (isoforms, contaminants, decoys) |

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
