---
name: bio-proteomics-spectral-libraries
category: Data Analysis
description: Builds and manages DIA spectral libraries as peptide query parameters (precursor m/z, a few fragment m/z plus relative intensities, normalized RT, optional CCS), covering experimental DDA, chromatogram, and in-silico predicted libraries via Koina-served Prosit, AlphaPeptDeep, MS2PIP, and DeepLC, with iRT/CiRT RT calibration, NCE tuning, format conversion (DIA-NN tsv/speclib/parquet, OpenSWATH pqp/TraML, Spectronaut, blib/dlib/elib), and library QC/merge. Use when generating, calibrating, converting, or merging a spectral library to drive a DIA search. Running the actual DIA search is dia-analysis; building from DDA identifications depends on peptide-identification; modified-peptide libraries route to ptm-analysis; quantifying the result is quantification.
tool_type: mixed
primary_tool: encyclopedia
license: MIT
author: GPTomics
---

## Version Compatibility

Reference examples tested with: koinapy 0.0.5+ (checked on 0.0.11), ms2pip 4.0+ (checked on 4.2.0),
deeplc 4.1+ (checked on 4.5.0), pandas 2.2+

Install: `pip install koinapy ms2pip deeplc pandas numpy scipy pyteomics`. CLI tools: EncyclopeDIA
(Java), EasyPQP/FragPipe for DDA libraries, OpenMS for `TargetedFileConverter` and
`OpenSwathDecoyGenerator` (checked on OpenMS 3.5.0).

deeplc dropped the `DeepLC()` class in 4.1 for a module-level API built on `psm_utils.PSM`/
`PSMList` (verified running on deeplc 4.5.0):

```python
from psm_utils import PSM, PSMList
import deeplc

cal_psms = PSMList(psm_list=[PSM(peptidoform=seq, spectrum_id=str(i), retention_time=rt)
                              for i, (seq, rt) in enumerate(calibration_pairs)])  # observed RT required
pred_psms = PSMList(psm_list=[PSM(peptidoform=seq, spectrum_id=str(i))
                               for i, seq in enumerate(peptides_to_predict)])
calibrated_rt = deeplc.predict_and_calibrate(pred_psms, psm_list_reference=cal_psms)
```

`deeplc.predict()` alone returns RT on the model's internal training-gradient units, not the
run's real RT -- always pair it with `deeplc.calibrate()`/`predict_and_calibrate()` against
observed-RT anchors, same rule as Koina iRT below. Modifications use ProForma notation via
`Peptidoform` (e.g. `'LGGNEQVTR[Phospho]'`), not the old MS2PIP `location|name` string.

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# DIA Spectral Libraries -- Query Parameters That Are Only as Good as Their Empirical Calibration

**"Build a spectral library for my DIA search"** -> Assemble a table of peptide query parameters (precursor m/z, top fragment m/z plus relative intensities, normalized RT, optional CCS), then calibrate the predicted RT/CCS to the actual gradient/instrument -- because a DIA library is not whole spectra and an uncalibrated prediction extracts every peak group at the wrong time.
- Python: `koinapy.Koina(...).predict(df)` for Prosit/AlphaPeptDeep/MS2PIP/UniSpec fragment intensities and iRT served from Koina
- Python: `deeplc.predict_and_calibrate(psms, psm_list_reference=cal_psms)` for RT prediction of any modification (module-level API, see Version Compatibility)
- Python: `ms2pip.predict_batch(psms, model='HCD')` for local fragment-intensity prediction
- CLI: EncyclopeDIA for empirical chromatogram libraries; EasyPQP/FragPipe for DDA-based libraries

Scope: this skill owns library generation (experimental, chromatogram, predicted, empirically-corrected), RT/CCS/NCE calibration, format conversion, and library QC/merge. Running the DIA search against the library is dia-analysis. Generating the DDA peptide identifications a DDA library is built from depends on peptide-identification. Modified-peptide and PTM-resolved library design routes to ptm-analysis. Quantifying and rolling up the search output is quantification. OUT OF SCOPE: acquisition-window design (fixed/variable/staggered/diaPASEF) and demultiplexing -- those belong to dia-analysis.

## The Single Most Important Modern Insight -- A Predicted Library Is Only as Good as Its Empirical Calibration

1. **A DIA library is peptide QUERY PARAMETERS, not spectra.** Each entry is a precursor m/z, a handful of fragment m/z with RELATIVE intensities, a normalized RT, and optionally CCS -- the inputs to extract and score a co-eluting fragment-chromatogram peak group, not a lookup spectrum to match (Gillet 2012). The fragments and their relative intensities are the discriminating content; absolute intensity is irrelevant.

2. **Fragment-intensity prediction is robust; predicted RT and CCS are in ARBITRARY model units and MUST be calibrated.** HCD fragmentation is reproducible at matched collision energy, so predicted relative intensities transfer across instruments. But predicted iRT is an arbitrary scale and real RT depends on the exact column, gradient, temperature, and mobile phase. The catastrophic error: drop a predicted library straight into a search without anchoring its RT to observed RT (via iRT/CiRT spike-in peptides or a GPF-DIA empirical pass). Every peak group is then extracted at the wrong time, selectivity collapses, and IDs silently disappear with no error. The same holds for AlphaPeptDeep CCS against the timsTOF's measured 1/K0.

3. **NCE (normalized collision energy) must match the predictor's training or intensities mismatch the real spectra.** Intensity predictors are conditioned on NCE. Feed a value that differs from the instrument's effective NCE and predicted intensities diverge from reality, losing sensitivity with no error. Do not blindly use NCE=30 from a tutorial -- scan candidate NCE values, predict, and pick the one maximizing spectral contrast/correlation against a few real spectra.

## Tool Taxonomy

| Library type / tool | Citation | Mechanism / role | When |
|---------------------|----------|------------------|------|
| Experimental DDA (EasyPQP, FragPipe) | -- | Consensus spectra from DDA runs of the same/pooled sample | Deep DDA already in hand; gold-standard real intensities |
| Chromatogram library (EncyclopeDIA, GPF) | Searle 2018 | GPF-DIA of pooled sample, narrow staggered windows -> empirical RT + real fragments in the actual LC | One project, maximum depth without fractionated DDA |
| In-silico predicted (Prosit, AlphaPeptDeep, MS2PIP+DeepLC) | Gessulat 2019; Zeng 2022 | Deep learning predicts fragment intensities + RT (+CCS) for the whole FASTA digest | No wet-lab library; the default modern route |
| Empirically-corrected predicted (EncyclopeDIA) | Searle 2020 | Predict whole-proteome library, search one GPF-DIA pass, rewrite intensities + RT with observed values | Non-model organisms, variant DBs; best of predicted + empirical |
| Prosit (intensity + iRT) | Gessulat 2019 | HCD/CID intensity conditioned on NCE; iRT model | Served via Koina; broad default predictor |
| AlphaPeptDeep (intensity + RT + CCS) | Zeng 2022 | Modular, retrainable; b/y plus mod neutral losses; predicts CCS | Full predicted library including ion mobility |
| MS2PIP (intensity) | -- | Fast HCD/CID/TMT/immuno intensity models; RT via DeepLC | Local prediction without a server; pairs with DeepLC |
| DeepLC (RT, any modification) | Bouwmeester 2021 | RT prediction for novel/modified peptides; needs calibration peptides | RT for peptidoforms carrying unseen modifications |
| UniSpec (NIST) | -- | Full-range intensity including internal/immonium ions | Available on Koina when richer fragment sets are needed |
| SpectraST (TPP) | -- | Legacy DDA consensus library builder | Legacy only; prefer EasyPQP/FragPipe instead |
| Acquisition window design / demux | -- | Fixed/variable/staggered/diaPASEF schemes | route OUT -> dia-analysis |

## Decision Tree by Scenario

| Scenario | Recommended | Why |
|----------|-------------|-----|
| No prior DDA, model organism, standard mods | Predicted library (Prosit or AlphaPeptDeep via Koina) + iRT calibration | Whole-proteome coverage, bounded search; calibrate RT to the gradient. Koina how-to: `references/koina-prediction.md` |
| Need ion mobility (timsTOF/diaPASEF) | AlphaPeptDeep (intensity + RT + CCS) | Only predictor here that emits CCS; calibrate CCS to measured 1/K0 |
| Non-model organism or custom/variant DB | Empirically-corrected predicted (Searle 2020) | One GPF-DIA pass rewrites predicted intensities/RT with observed values |
| Maximum depth for one project, have pooled sample | Chromatogram library (EncyclopeDIA + GPF) | Empirical RT and real fragmentation in the actual LC |
| Deep fractionated DDA already acquired | Experimental DDA library (EasyPQP/FragPipe) | Real consensus spectra; gold-standard intensities |
| Library for OpenSWATH | Any source, then OpenSwathDecoyGenerator | OpenSWATH needs decoys IN the library; target-only has no null. Conversion and decoy requirements: `references/format-conversion.md` |
| Modified/PTM peptidoforms required | Include mods in digest; DeepLC for RT of unseen mods | route to ptm-analysis for PTM-resolved design |

Default when uncertain: a Koina-served predicted library (Prosit intensity + iRT) with explicit iRT/CiRT RT calibration and an NCE scan, exported to the search engine's native format.

### Calibrate iRT to Observed RT

**Goal:** Map arbitrary-unit predicted iRT onto the run's real RT so peak groups extract at the right time.

**Approach:** Spike or detect anchor peptides (11 Biognosys iRT peptides, or CiRT endogenous peptides when no spike-in exists), fit a regression from library iRT to observed RT, and require a tight fit before trusting it. A global linear fit fails on nonlinear gradients -- fall back to LOWESS.

```python
import numpy as np
from scipy import stats

IRT_PEPTIDES = {'LGGNEQVTR': -24.92, 'GAGSSEPVTGLDAK': 0.00, 'VEATFGVDESNAK': 12.39,
                'YILAGVENSK': 19.79, 'TPVISGGPYEYR': 28.71, 'TPVITGAPYEYR': 33.38,
                'DGLDAASYYAPVR': 42.26, 'ADVTPADFSEWSK': 54.62, 'GTFIIDPGGVIR': 70.52,
                'GTFIIDPAAVIR': 87.23, 'LFLQFGAQGSPFLK': 100.00}

R2_MIN = 0.95  # below this the RT alignment is untrustworthy and extraction windows are misplaced

def fit_irt_to_rt(anchor_irt, observed_rt):
    slope, intercept, r, _, _ = stats.linregress(anchor_irt, observed_rt)
    if r ** 2 < R2_MIN:
        raise ValueError(f'iRT fit R^2={r**2:.3f} < {R2_MIN}; gradient may be nonlinear, use LOWESS')
    return lambda irt: slope * irt + intercept
```

## Reference Files

| File | Read when |
|------|-----------|
| `references/koina-prediction.md` | Predicting fragment intensities and iRT for a peptide list via Koina (Prosit/AlphaPeptDeep/MS2PIP/UniSpec), including peptide validation before the request |
| `references/format-conversion.md` | Converting between DIA-NN, OpenSWATH and Spectronaut; building an OpenSWATH transition TSV (`scripts/build_openswath_tsv.py`, `scripts/spectronaut_to_diann.py`); OpenSwathDecoyGenerator requirements and decoy method choice |
| `references/library-qc-merge.md` | Reporting library size and merging libraries on the full transition key |

## Per-Method Failure Modes

### Predicted RT/CCS used without calibration
**Trigger:** A predicted library is searched directly, RT column straight from the model.
**Mechanism:** Predicted iRT/CCS are arbitrary model units; real RT depends on column/gradient/temperature.
**Symptom:** Drastic ID loss with no error; peak groups extracted at the wrong time.
**Fix:** Fit iRT/CiRT anchors (R^2 > 0.95) or run a GPF-DIA empirical correction (Searle 2020) before searching.

### NCE mismatch
**Trigger:** A fixed collision energy (often 30) reused across instruments/methods.
**Mechanism:** Intensity predictors are NCE-conditioned; wrong NCE shifts predicted relative intensities.
**Symptom:** Quiet sensitivity loss; fewer confident peak groups than expected.
**Fix:** Scan candidate NCE values, predict, and pick the one maximizing spectral contrast against real spectra.

### Missing decoys for OpenSWATH
**Trigger:** A target-only library handed to OpenSWATH.
**Mechanism:** Peptide-centric scoring needs a decoy null; OpenSWATH does not invent one.
**Symptom:** FDR cannot be estimated or is meaningless.
**Fix:** Run OpenSwathDecoyGenerator to append decoys; do NOT also supply decoys to DIA-NN/Spectronaut, which generate their own.

### OpenSwathDecoyGenerator silently generates 0 decoys or merges peptides
**Trigger:** Transition list has placeholder/approximate `ProductMz`, lacks a literal `Annotation` column, or lacks `transition_group_id`, yet converts cleanly via TargetedFileConverter.
**Mechanism:** The decoy algorithm matches target and decoy fragments by annotation and real m/z, so without both it pairs no fragment and drops every candidate. Without a grouping column TargetedFileConverter merges peptides sharing a charge into one group.
**Symptom:** "Number of decoy peptides: 0" and a hard failure at the 80% threshold check; or, with no error at all, fewer `<Peptide>` groups in the TraML than precursors in the TSV.
**Fix:** Supply `Annotation`, real theoretical fragment m/z and `transition_group_id`, and compare the `<Peptide>` count to the precursor count after conversion -- see `references/format-conversion.md`.

### Modification mismatch between library and data
**Trigger:** Library lacks the sample's variable mods, or carries too many.
**Mechanism:** A library without phospho/ox cannot find those peptidoforms; too many variable mods explode the search space and inflate FDR.
**Symptom:** Missing modified peptides, or inflated IDs.
**Fix:** Match library modifications to the biology; route PTM-resolved design to ptm-analysis.

### Naive merge dropping transitions
**Trigger:** Dedup keyed on (sequence, fragment-type, fragment-number) only.
**Mechanism:** Distinct transitions can share those three fields but differ in precursor or fragment charge.
**Symptom:** Quietly thinner transition lists; weaker peak-group scoring.
**Fix:** Key dedup on (modified-sequence, precursor-charge, fragment-type, fragment-number, fragment-charge).

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| 6 fragments per precursor | OpenSWATH/EncyclopeDIA defaults | Enough for confident peak-group scoring; more invites interference |
| Fragment m/z > precursor m/z, and > ~200 | Practice | Avoids the low-mass region dense with shared/uninformative ions |
| Library FDR 1% (peptide and protein) | EasyPQP defaults | A dirty library poisons every downstream search |
| iRT regression R^2 > 0.95 | Practice | Below this RT alignment is untrustworthy and windows misplace |
| NCE chosen by spectral-contrast scan | Practice | Matches the predictor's training to the instrument's effective NCE |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| ConnectionError on proteomicsdb.org/prosit/api/predict | The old Prosit endpoint is dead | Use Koina: `from koinapy import Koina; Koina('Prosit_2019_intensity', 'koina.wilhelmlab.org:443')` |
| ImportError: cannot import name Predictor from ms2pip | No Predictor class in ms2pip v4 | Call module-level `ms2pip.predict_batch(psms, model='HCD')` returning ProcessingResult objects |
| koinapy TypeError on constructor/columns | Constructor signature and column names vary by version | Verify with `help(Koina)`; inputs are typically `peptide_sequences`, `precursor_charges`, `collision_energies` |
| DeepLC RT all near constant or on the wrong scale | `deeplc.predict()` called without calibration | Use `deeplc.calibrate()` + `deeplc.predict()`, or `deeplc.predict_and_calibrate(psms, psm_list_reference=cal_psms)`; mods as ProForma via `Peptidoform`, not MS2PIP `location\|name` |
| Extraction at wrong time, ID collapse | Predicted RT not calibrated to the gradient | Fit iRT/CiRT anchors or run GPF-DIA empirical correction before searching |
| OpenSWATH FDR meaningless | Target-only library, no decoys | Append decoys with OpenSwathDecoyGenerator |
| OpenSwathDecoyGenerator: "Number of decoy peptides: 0" / below 80% threshold, or fewer `<Peptide>` groups than precursors | Missing `Annotation`, placeholder `ProductMz`, or no `transition_group_id` | See the OpenSwathDecoyGenerator failure mode above |
| Decoy peptide sequences differ between identical `-method shuffle` runs | Shuffle decoy generation has no seed flag; this is expected, not a bug | Use `-method pseudo-reverse` for reproducible decoys |
| ms2pip.predict_batch appears to hang on first call ("Model hash not recognized." then nothing) | The default `model='HCD'` (= HCD2021) downloads two XGBoost files to `~/.ms2pip` on first use -- 66MB + 847MB, confirmed via `Content-Length` -- with no progress output, no timeout, and no resume (a killed download restarts from 0, not where it left off) | Let it finish once with network access (~10+ min on a slow link), pre-populate `model_dir` from a machine that already has it cached, or pass a smaller model (`model='HCD2019'`, ~17MB total, confirmed working end-to-end in ~45s) if HCD2021's extra accuracy isn't needed |
| Fewer transitions than expected after merge | Dedup key missed charges | Key on the full five-field transition key (`references/library-qc-merge.md`) |

## References

- Gillet LC, Navarro P, Tate S, et al. Targeted data extraction of the MS/MS spectra generated by data-independent acquisition: a new concept for consistent and accurate proteome analysis. *Mol Cell Proteomics* 2012;11(6):O111.016717.
- Searle BC, Pino LK, Egertson JD, et al. Chromatogram libraries improve peptide detection and quantification by data independent acquisition mass spectrometry. *Nat Commun* 2018;9:5128.
- Gessulat S, Schmidt T, Zolg DP, et al. Prosit: proteome-wide prediction of peptide tandem mass spectra by deep learning. *Nat Methods* 2019;16(6):509-518.
- Searle BC, Swearingen KE, Barnes CA, et al. Generating high quality libraries for DIA MS with empirically corrected peptide predictions. *Nat Commun* 2020;11:1548.
- Bouwmeester R, Gabriels R, Hulstaert N, Martens L, Degroeve S. DeepLC can predict retention times for peptides that carry as-yet unseen modifications. *Nat Methods* 2021;18(11):1363-1369.
- Zeng WF, Zhou XX, Willems S, et al. AlphaPeptDeep: a modular deep learning framework to predict peptide properties for proteomics. *Nat Commun* 2022;13:7238.

## Related Skills

- dia-analysis - Run the DIA search against the library and choose the q-value context
- peptide-identification - Generate the DDA identifications a DDA library is built from
- ptm-analysis - Design PTM-resolved and modified-peptide libraries
- quantification - Summarize and roll up the search output to protein abundances
