## Per-Tool Failure Modes

### ADMETlab 3.0 -- out-of-distribution prediction

**Trigger:** Chemistry materially unlike the service's documented training/applicability domain, such as metal-containing complexes, many peptides, PROTACs, or unusual macrocycles.

**Mechanism:** ADMETlab training set is drug-like organic molecules. Predictions on PROTACs, macrocycles, peptides extrapolate.

**Symptom:** High reported uncertainty (where the hosted response supplies it), disagreement with neighbors or orthogonal models, or unstable conclusions under reasonable preprocessing.

**Fix:** For an ADMETlab response, inspect its modeled-endpoint uncertainty using the current service contract. For the offline ADMET-AI route, run `assess_applicability_domain` against a relevant local reference set before prediction; it rejects metals/inorganics and flags low-similarity or oversized organic chemistry for manual review. For PROTACs / macrocycles, prefer literature-derived experimental data.

### hERG D-MPNN -- training data bias

**Trigger:** Compound is novel chemotype not in training set (drug-like but in unexplored region).

**Mechanism:** D-MPNN learns local chemical features; for genuinely new scaffolds, extrapolation is unreliable.

**Symptom:** Model predicts hERG- (false negative) for compound that experimentally inhibits.

**Fix:** Use ensemble + applicability-domain assessment (kNN distance, ensemble variance). If kNN distance to training set > P95, treat prediction as low-confidence. For the bundled offline route, use `fit_similarity_threshold` and `assess_applicability_domain` as a transparent reference-similarity screen; it is not a substitute for model-specific AD validation.

### CYP3A4 inhibitor + substrate ambiguity

**Trigger:** Model trained on either inhibitor OR substrate; predictions confused.

**Mechanism:** CYP3A4 inhibitors and substrates have similar SAR; many compounds are both.

**Symptom:** Both classes report > 0.5.

**Fix:** Two separate models (inhibitor model, substrate model); compounds that score high in both are flagged for in vitro confirmation.

### SwissADME -- no API

**Trigger:** Wanting to batch programmatically.

**Mechanism:** SwissADME's terms restrict automated crawler/data-retrieval access, and no public API is documented.

**Symptom:** No programmatic access; manual web upload only.

**Fix:** Use a currently documented programmatic service and follow its access policy; for ADMETlab 3.0, verify the live API tutorial before writing a client.

### PAINS as a kill filter

**Trigger:** Treating PAINS_A match as a categorical exclusion.

**Mechanism:** PAINS is calibrated against HTS assay-interference; matches do NOT predict failed drug development.

**Symptom:** Library purged of valid leads (curcumin analogs, polyphenol natural products).

**Fix:** Flag PAINS for orthogonal-assay confirmation; do not exclude pre-emptively. See substructure-search for details.

### Class-imbalanced AMES dataset

**Trigger:** Training/predicting AMES mutagenicity.

**Mechanism:** Public AMES datasets can be imbalanced and differ in assay definition and curation; aggregate accuracy can therefore be misleading.

**Symptom:** Model reports high accuracy but predicts negative for all.

**Fix:** Report class balance and use suitable metrics such as PR-AUC, ROC-AUC, MCC, or balanced accuracy. Compare class weighting or resampling inside training folds without leaking validation/test data.

## Reconciliation Across Models

When ADMETlab, ProTox-3.0, and an independently trained chemprop model disagree on hERG:
- All predict hERG+ -> higher concern; plan in vitro patch-clamp
- Results disagree -> inspect applicability domains, activity thresholds, and assay definitions before deciding
- All predict hERG- -> lower concern, but still consider in vitro screening for clinical candidates and novel chemotypes
- Do not count correlated models as independent evidence merely because they are hosted by different services

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| ADMETlab API timeout | Service load, payload, or current quota | Follow live batch limits; retry with backoff and record failures |
| chemprop training overfits | Random split | Use scaffold split (`--split scaffold_balanced`) |
| hERG prediction 50/50 | Out-of-distribution | Check applicability domain |
| QED calculation fails | Molecule is missing, unsanitized, or unsupported | Reject parse failures; sanitize inputs and handle calculation exceptions |
| ProTox endpoints missing | Web scrape uses CSS selector | Use formal API |
| BBB+ true but TPSA > 90 | Different BBB model | Use the simple physicochemical screen as an orthogonal heuristic |
| Predictions inconsistent across runs | Data split and PyTorch randomness for chemprop ensemble | Use `--data-seed 42 --pytorch-seed 42` and reuse the saved model; verify flags against the installed Chemprop release |
| Calibration mismatch | DL native probabilities not calibrated | Apply Platt scaling on validation set |
