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
