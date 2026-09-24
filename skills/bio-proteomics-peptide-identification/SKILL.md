---
name: bio-proteomics-peptide-identification
description: Peptide-spectrum matching from MS/MS with target-decoy FDR control, framing identification confidence as a property of a ranked list (q-value/PEP) rather than a raw engine score (XCorr or SpecEValue). Covers the verified sequence-database search engines Comet, MS-GF+ and Sage; concatenated vs separate target-decoy competition; PEP vs q-value; the multi-level FDR cascade; Percolator rescoring; and pyOpenMS SimpleSearchEngineAlgorithm + FalseDiscoveryRate. Use when identifying peptides from tandem mass spectra and deciding what FDR threshold to act on. Protein grouping and protein-level FDR are protein-inference; PTM site localization and open-search follow-up are ptm-analysis; DIA peptide-centric scoring is dia-analysis; intensity quant is quantification.
tool_type: mixed
primary_tool: pyOpenMS
license: MIT
author: GPTomics
---

## Version Compatibility

Reference examples tested with: pyOpenMS 3.5.0, pandas 2.2+, numpy 1.26+ (pyOpenMS 3.5 takes a `PeptideIdentificationList`, not a plain Python list, for peptide IDs). Command-line route checked on Sage 0.14.6, Comet 2026.02 rev.2, MS-GF+ v2024.03.26, Percolator 3.09.0, OpenMS 3.5.0 `DecoyDatabase`, Java 17.

Install: `pip install pyopenms pandas numpy`. Command-line tools: `sage`, `comet`, `MSGFPlus` (`java -jar`), `percolator`, OpenMS `DecoyDatabase`; vendor raw -> mzML with `msconvert` (ProteoWizard) or ThermoRawFileParser.

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Peptide Identification -- Confidence Is a Property of a Ranked List, Not a Single PSM

**"Identify peptides from my MS/MS spectra"** -> Match tandem mass spectra against a protein database, then control false discovery rate by target-decoy competition and act on a q-value -- because a raw match score is meaningless in isolation; only the list-level error rate is interpretable.
- Python: `pyopenms.SimpleSearchEngineAlgorithm().search(...)` for in-process database search, `FalseDiscoveryRate` for q-values
- CLI: `comet`, `sage`, `MSGFPlus` for high-throughput database searching; `percolator` for rescoring

Scope: this skill owns spectrum-to-peptide matching and PSM/peptide-level FDR. Protein grouping and protein-level (picked) FDR -> protein-inference. PTM site localization and open-search mod discovery follow-up -> ptm-analysis. DIA peptide-centric extraction and scoring -> dia-analysis. FDR-filtered IDs feeding intensities -> quantification. mzML/raw loading -> data-import. OUT OF SCOPE: protein inference, PTM localization scoring, DIA peptide-centric pipelines, label-free/TMT quantification.

## The Single Most Important Modern Insight -- A q-value Is a Verdict on the List, a Raw Score Is Not Even Comparable

1. **Identification confidence is a property of a ranked LIST controlled by target-decoy competition, never a property of one PSM.** The number to act on is a q-value (list-level) or PEP (per-PSM), NOT the engine's raw score. XCorr (Comet) and SpecEValue (MS-GF+) live on different scales, are charge- and length-dependent, and are frequently not even monotone in true probability within a single engine -- which is why rescoring with Percolator exists. "1% FDR" answers "what fraction of the list I keep is wrong," NOT "I am 99% sure of this one ID." The catastrophic error is thresholding on a raw score, or comparing scores across engines.

2. **A q-value is valid only if (a) the decoy DB is a faithful null, (b) targets and decoys competed in ONE concatenated search, and (c) there are enough PSMs for the decoy count to be stable.** Generate decoys at the PROTEIN level then digest (so decoy peptides obey the same enzyme rules), matching the target in size and composition. Concatenated competition (one best hit per spectrum) gives FDR = (#decoys above threshold + 1) / (#targets above threshold) -- one decoy above threshold estimates one false target, and the +1 keeps small lists honest. Elias & Gygi's 2 * #decoy / (#target + #decoy) is the older, conservative form of the same concatenated-search estimate, counted over the whole target+decoy list. Separate target/decoy searches (no per-spectrum competition) instead need pi0 * #decoy / #target (Kall, Storey, MacCoss & Noble 2008) or the refined mix-max estimator (Keich, Kertesz-Farkas & Noble 2015). Applying 2d/(t+d) to separate searches over-estimates FDR and throws away identifications (synthetic test: 2.06% estimated vs 0.62% true).

3. **PEP and q-value answer different questions; filtering at "PEP <= 0.01" is far stricter than "q <= 0.01."** PEP (posterior error probability, local FDR) is the probability that THIS PSM is wrong; q-value is the FDR of the list cut at this PSM. FDR is the average of PEP over the accepted set (Kall 2008). The worst PSM in a 1%-FDR list typically has a PEP of 10-50% (`examples/fdr_filtering.py` prints the q-value and PEP cuts side by side on a demo table). Use q-value for list cutoffs; use PEP only for per-ID decisions (e.g. picking one PTM site). And PSM-FDR at 1% does NOT give 1% peptide-FDR or 1% protein-FDR -- each level needs its own estimation; hand protein-level control to protein-inference.

## The FDR Vocabulary, Precisely

- **FDR**: the expected proportion of false positives among ALL accepted items at a threshold -- a property of the whole list.
- **q-value**: the minimum FDR at which a given PSM is still accepted; monotone after taking the running minimum from the bottom of the ranked list. Filter on q <= 0.01.
- **PEP (local FDR)**: the probability that THIS PSM is wrong given its score. Local, per-PSM; FDR is the integral of PEP over the accepted set (Kall 2008, "two sides of the same coin").
- **The estimator must match the search mode.** Concatenated target-decoy competition (TDC): FDR = (#decoy + 1) / #target (one best hit per spectrum already resolves the competition; Elias-Gygi's 2 * #decoy / (#target + #decoy) is the older, conservative whole-list form for this same composite search). Separate target and decoy searches: pi0 * #decoy / #target (Kall et al. 2008; pi0 = 1 is valid but conservative), or the refined mix-max estimator (Keich, Kertesz-Farkas & Noble 2015). Mix-max is a distinct, calibrated-score procedure for the separate-search setting.

## Tool Taxonomy

| Tool / method | Citation | Mechanism / role | When |
|---|---|---|---|
| Comet | Eng 2013 | XCorr + E-value; SEQUEST lineage, open-source | Robust default, TPP pipelines; pairs with Percolator |
| MS-GF+ | Kim & Pevzner 2014 | SpecEValue via generating-function DP | Calibrated cross-instrument E-value; ETD/CID, low-res, non-standard enzymes |
| Sage | Lazear 2023 | hyperscore-style, Rust, rescoring-native | Modern scalable open-source pipelines; emits Percolator-ready features |
| Percolator | Kall 2007 | semi-supervised SVM re-rank on decoy negatives | Boost IDs at fixed FDR; non-tryptic/PTM/large search spaces |
| Spectral-library search | -- | match empirical reference spectra (intensity + RT) | Faster/more specific for known peptides -> spectral-libraries |
| Protein grouping / protein FDR | Savitski 2015; The 2022 | picked / picked-group FDR | route OUT -> protein-inference |
| PTM site localization | -- | per-site PEP, localization scoring | route OUT -> ptm-analysis |

## Decision Tree by Scenario

| Scenario | Recommended | Why |
|---|---|---|
| Standard DDA, clean FDR, scriptable | Comet or Sage + Percolator at q <= 0.01 (`references/cli_route.md`) | well-validated; rescoring boosts IDs at fixed FDR |
| Search or filter inside Python, no external engine | pyOpenMS `SimpleSearchEngineAlgorithm` + `FalseDiscoveryRate` (`references/pyopenms.md`) | in-process; output already target/decoy annotated |
| Have a PSM table from any engine, need q-values | concatenated table, or pi0 * D / T for separate searches (`references/fdr_from_tables.md`) | the estimator must match the search mode |
| Cross-instrument / varied fragmentation / odd enzyme | MS-GF+ | SpecEValue is calibrated so a threshold means the same everywhere |
| Discover unknown PTMs / mass shifts | ptm-analysis | open-search setup and localization follow-up belong there |
| Huge dataset, reproducible, cloud-scale | Sage (rescoring-native) | Rust speed; emits Percolator features directly |
| All-in-one search and quantification | quantification | this Skill ends at FDR-filtered IDs; intensity analysis belongs there |
| Non-tryptic (immunopeptidomics, degradomics) | a verified engine + Percolator | rescoring gains are largest where search space explodes |
| Few PSMs (single-protein pulldown) | do NOT trust decoy FDR; inspect spectra manually | decoy counts too noisy below ~hundreds of PSMs |
| Need per-site / per-ID confidence | act on PEP, not q-value | q-value is list-level; PEP is local |

Default when uncertain: concatenated target-decoy search with Comet or Sage, rescore with Percolator, filter at q <= 0.01, and hand protein-level FDR to protein-inference.

### Build the Concatenated Target-Decoy Database

**Goal:** Turn a target-only FASTA into the concatenated target+decoy database that every FDR number below depends on.

**Approach:** Reverse at the PROTEIN level with the peptide N/C termini held fixed, so decoy peptides obey the same enzyme rules and match the targets in length and amino-acid composition, at 1:1. Sage generates decoys internally and takes the TARGET-ONLY FASTA; Comet and MS-GF+ need the concatenated file. Include contaminants in the input before reversing, so contaminant decoys exist too. **The decoy tag is where this Skill silently fails:** every tool defaults to a different string, and a mismatch means zero decoys are found and every PSM is reported "at 1% FDR" with no error.

```bash
# OpenMS 3.5.0; -method shuffle is the alternative when reversal makes
# decoy peptides that are palindromes of real ones
DecoyDatabase -in human_plus_contaminants.fasta -out human_target_decoy.fasta \
  -decoy_string DECOY_ -decoy_string_position prefix \
  -method reverse -enzyme Trypsin -threads 8

grep -c '^>' human_target_decoy.fasta        # must be 2x the target count
grep -c '^>DECOY_' human_target_decoy.fasta  # must equal the target count
```

| Tool | Default decoy tag | Set it with |
|---|---|---|
| OpenMS `DecoyDatabase` / `PeptideIndexing` | `DECOY_` | `-decoy_string`, `decoy_string` |
| Sage | `rev_` (lower case) | `decoy_tag` in the JSON, with `generate_decoys: true` |
| Comet | `DECOY_` | `decoy_prefix`; `decoy_search = 0` when the DB already holds them |
| MS-GF+ | `XXX_` | `-decoy`; `-tda 0` when the DB already holds them |
| Percolator | reads the pin `Label` column | `-P` only for `--picked-protein` |

Checked on OpenMS 3.5.0: 31,437 UniProt entries (human + yeast + E. coli + contaminants) became 62,874 with 31,437 `DECOY_` entries in 5 s.

## Reference Files

Read the file that matches the task; `SKILL.md` keeps only what every request needs.

| File | Read when |
|---|---|
| `references/cli_route.md` | running Sage, Comet or MS-GF+ from the command line; rescoring the `.pin` with Percolator; pooling runs (`examples/dda_search.sh` runs it end to end) |
| `references/pyopenms.md` | searching or FDR-filtering inside Python with `SimpleSearchEngineAlgorithm`, `PeptideIndexing`, `FalseDiscoveryRate` |
| `references/fdr_from_tables.md` | q-values from a PSM table of any engine, or from separate target and decoy searches (`examples/separate_search_fdr.py`) |
| `references/citations.md` | the full citation list |

Runnable code: `scripts/pyopenms_search.py`, `scripts/pyopenms_fdr.py`, `scripts/table_fdr.py` (each with `--help`); end-to-end examples in `examples/`.

## Per-Method Failure Modes

### Concatenated vs separate FDR formula mismatch
**Trigger:** running the concatenated-search code (`scripts/table_fdr.py`) on a merged table from separate searches without per-spectrum competition, or applying Elias-Gygi's 2*decoy/(target+decoy) to separate searches.
**Mechanism:** Elias-Gygi's factor 2 counts decoys in the combined target+decoy list of a concatenated search; in separate searches every spectrum gets both a target and a decoy hit, so the decoy count estimates false targets directly, scaled by pi0 (the fraction of target PSMs that are incorrect).
**Symptom:** mis-estimated FDR; 2d/(t+d) on separate searches over-estimates it (synthetic test: 2.06% vs 0.62% true, about a quarter of IDs lost).
**Fix:** confirm the search mode; concatenated -> (#decoy + 1)/#target; separate -> pi0 * #decoy/#target (Kall et al. 2008; code in `references/fdr_from_tables.md` and `examples/separate_search_fdr.py`) or the mix-max estimator (Keich, Kertesz-Farkas & Noble 2015). In Percolator, mix-max is the default for separate-search input and `-Y`/`--post-processing-tdc` selects target-decoy competition instead; concatenated input forces TDC automatically.

### Thresholding on raw engine score
**Trigger:** filtering on XCorr or another raw engine score, or comparing scores from two engines.
**Mechanism:** scores are uncalibrated, charge/length-dependent, and not monotone in true probability.
**Symptom:** different cutoffs admit different real FDRs; cross-engine merges nonsensical.
**Fix:** always convert to q-value (or SpecEValue/PEP) first; rescore with Percolator.

### Decoy FDR on too few PSMs
**Trigger:** reporting "0% FDR" from a single-protein pulldown or tiny PSM list.
**Mechanism:** the decoy count is a noisy Poisson-like estimate; zero observed decoys does not mean zero false targets.
**Symptom:** spuriously confident IDs from small experiments.
**Fix:** below ~hundreds of PSMs, inspect spectra manually; do not act on the decoy q-value.

### Open-search results used for clean FDR or quant
**Trigger:** taking IDs from a wide-window (-150..+500 Da) search as final, FDR-controlled results.
**Mechanism:** wide windows admit "free" mass shifts that inflate random matches; the target-decoy null differs per mass-shift bin.
**Symptom:** inflated, unreliable FDR on open-search output.
**Fix:** treat open search as discovery; follow with a closed search restricted to the discovered mods -> ptm-analysis.

### Rescoring overfitting
**Trigger:** custom features that leak label information, or training without proper cross-validation.
**Mechanism:** the model learns the decoys, making rescored FDR optimistic.
**Symptom:** ID counts jump but downstream validation fails.
**Fix:** use Percolator's default cross-validation; validate with entrapment for high-stakes claims (Wen 2025).

### DIA tool FDR taken at face value
**Trigger:** trusting a DIA tool's reported 1% peptide/protein FDR.
**Mechanism:** entrapment shows several DIA tools do not reliably control FDR (Wen 2025).
**Symptom:** real error rate exceeds the reported FDR.
**Fix:** validate with entrapment for high-stakes DIA claims -> dia-analysis.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|---|---|---|
| Precursor tolerance 10-20 ppm (high-res Orbitrap) | -- | matches FT mass accuracy; tighter = fewer random candidates at fixed FDR |
| Fragment tolerance 0.02 Da (HCD Orbitrap) / 0.6 Da (ion-trap CID) | -- | instrument-dependent; 0.6 Da on Orbitrap discards resolving power |
| Missed cleavages 2 | -- | covers incomplete trypsin digestion without exploding search space |
| PSM/peptide FDR 1% (q <= 0.01) | Elias & Gygi 2007 | community standard; list-level error, not per-PSM |
| Decoy:target ratio 1:1 | Elias & Gygi 2007 | standard; unequal ratios need formula correction |
| Min PSMs for trustworthy decoy FDR: hundreds+ | -- | below this the decoy count is too noisy |
| Variable mods per peptide <= 2-3 | -- | each variable mod multiplies search space and random-match rate |

## Common Errors

| Error / symptom | Cause | Solution |
|---|---|---|
| pyOpenMS "search" returns peptides but never scores spectra | used `ProteaseDigestion`, which only digests a FASTA | use `SimpleSearchEngineAlgorithm().search(mzML, fasta, protein_ids, peptide_ids)` |
| `TypeError: Argument 'pep_ids' has incorrect type (expected ...PeptideIdentificationList, got list)` or `can not handle type` | pyOpenMS 3.5+ needs a `PeptideIdentificationList` | `peptide_ids = PeptideIdentificationList()`; protein_ids FIRST: `IdXMLFile().load(path, protein_ids, peptide_ids)` |
| `RuntimeError: Meta value 'target_decoy' does not exist` from `FalseDiscoveryRate` | decoys not annotated (e.g. idXML from another engine) | run `PeptideIndexing` with matching `decoy_string` first |
| Every PSM passes 1% FDR, or `ValueError: no decoy PSMs recognised` | decoy prefix not matched (Sage writes lowercase `rev_`) | compare prefixes lower-cased; check `psms['protein'].str[:6].value_counts()` |
| Empty 1% list from `scripts/table_fdr.py` | lower-is-better score (E-value, SpecEValue) used as `score` | use `-log10(E-value)` |
| All q-values 0 from a hand-rolled table | no +1 correction on a list with zero decoys | use (decoys + 1)/targets; a tiny list cannot reach 1% |
| Percolator q-method mismatched to search mode | mix-max is the default for separate-search input | for separate searches, mix-max (default) or `-Y`/`--post-processing-tdc` for target-decoy competition; concatenated input forces TDC automatically; use `--picked-protein` for protein FDR |
| Percolator's 1% list is far too big or too small when cut by column index | Percolator writes a `filename` column only when the pin has one (Sage yes, Comet no), shifting `q-value` between columns 3 and 4 | locate `q-value` by header name, never by a fixed index |
| Decoy count is twice the target count after a search | `-tda 1` (MS-GF+) or `generate_decoys: true` (Sage) run against a database that already contains decoys | use `-tda 0` / `generate_decoys: false` with a concatenated DB, or feed the target-only FASTA and let the engine make them |
| 1% PSM FDR assumed to give 1% protein FDR | each level needs its own estimation | estimate protein-level (picked) FDR -> protein-inference |
| "PEP <= 0.01" returns far fewer IDs than expected | PEP is per-PSM and far stricter than q-value | filter list cutoffs on q-value; reserve PEP for per-ID decisions |

## Related Skills

- protein-inference - Group peptides to protein groups and control protein-level (picked) FDR
- ptm-analysis - Open/variable-mod search follow-up and per-site PTM localization
- dia-analysis - DIA peptide-centric extraction and scoring; entrapment FDR validation
- quantification - FDR-filtered IDs feed label-free/TMT intensity quantification
- spectral-libraries - Empirical and predicted spectral-library search as an ID alternative
- data-import - Load mzML/raw MS data before identification
- database-access/uniprot-access - Build the target FASTA (canonical vs isoform, contaminants)
