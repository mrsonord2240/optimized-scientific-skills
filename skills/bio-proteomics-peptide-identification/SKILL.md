---
name: bio-proteomics-peptide-identification
description: Peptide-spectrum matching from MS/MS with target-decoy FDR control, framing identification confidence as a property of a ranked list (q-value/PEP) rather than a raw engine score (XCorr, hyperscore, Andromeda, SpecEValue). Covers sequence-database search engines (Comet, MS-GF+, MSFragger, Sage, MaxQuant, MetaMorpheus), concatenated vs separate target-decoy competition, PEP vs q-value, the multi-level FDR cascade, open/mass-tolerant search, rescoring (Percolator, mokapot, MS2Rescore), and pyOpenMS SimpleSearchEngineAlgorithm + FalseDiscoveryRate. Use when identifying peptides from tandem mass spectra and deciding what FDR threshold to act on. Protein grouping and protein-level FDR are protein-inference; PTM site localization is ptm-analysis; DIA peptide-centric scoring is dia-analysis; intensity quant is quantification.
tool_type: mixed
primary_tool: pyOpenMS
license: MIT
author: GPTomics
---

## Version Compatibility

Reference examples tested with: pyOpenMS 3.5.0, pandas 2.2+, numpy 1.26+ (pyOpenMS 3.5 takes a `PeptideIdentificationList`, not a plain Python list, for peptide IDs). Command-line route checked on Sage 0.14.6, Comet 2026.02 rev.2, MS-GF+ v2024.03.26, Percolator 3.09.0, OpenMS 3.5.0 `DecoyDatabase`, Java 17.

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Peptide Identification -- Confidence Is a Property of a Ranked List, Not a Single PSM

**"Identify peptides from my MS/MS spectra"** -> Match tandem mass spectra against a protein database, then control false discovery rate by target-decoy competition and act on a q-value -- because a raw match score is meaningless in isolation; only the list-level error rate is interpretable.
- Python: `pyopenms.SimpleSearchEngineAlgorithm().search(...)` for in-process database search, `FalseDiscoveryRate` for q-values
- CLI: `comet`, `msfragger`, `sage`, `MSGFPlus` for high-throughput database searching, `percolator`/`mokapot` for rescoring
- R: `mzID::mzID()` + `flatten()` or `mzR::openIDfile()` + `psms()` to read mzIdentML search results

Scope: this skill owns spectrum-to-peptide matching and PSM/peptide-level FDR. Protein grouping and protein-level (picked) FDR -> protein-inference. PTM site localization and open-search mod discovery follow-up -> ptm-analysis. DIA peptide-centric extraction and scoring -> dia-analysis. FDR-filtered IDs feeding intensities -> quantification. mzML/raw loading -> data-import. OUT OF SCOPE: protein inference, PTM localization scoring, DIA peptide-centric pipelines, label-free/TMT quantification.

## The Single Most Important Modern Insight -- A q-value Is a Verdict on the List, a Raw Score Is Not Even Comparable

1. **Identification confidence is a property of a ranked LIST controlled by target-decoy competition, never a property of one PSM.** The number to act on is a q-value (list-level) or PEP (per-PSM), NOT the engine's raw score. XCorr (Comet), hyperscore (MSFragger/X!Tandem), Andromeda score (MaxQuant), and SpecEValue (MS-GF+) live on different scales, are charge- and length-dependent, and are frequently not even monotone in true probability within a single engine -- which is exactly why rescoring (Percolator/mokapot) exists. "1% FDR" answers "what fraction of the list I keep is wrong," NOT "I am 99% sure of this one ID." The catastrophic error is thresholding on a raw score, or comparing scores across engines.

2. **A q-value is valid only if (a) the decoy DB is a faithful null, (b) targets and decoys competed in ONE concatenated search, and (c) there are enough PSMs for the decoy count to be stable.** Generate decoys at the PROTEIN level then digest (so decoy peptides obey the same enzyme rules), matching the target in size and composition. Concatenated competition (one best hit per spectrum) gives FDR = (#decoys above threshold + 1) / (#targets above threshold) -- one decoy above threshold estimates one false target, and the +1 keeps small lists honest. Elias & Gygi's 2 * #decoy / (#target + #decoy) is the older, conservative form of the same concatenated-search estimate, counted over the whole target+decoy list. Separate target/decoy searches (no per-spectrum competition) instead need pi0 * #decoy / #target (Kall, Storey, MacCoss & Noble 2008) or the refined mix-max estimator (Keich, Kertesz-Farkas & Noble 2015). Applying 2d/(t+d) to separate searches over-estimates FDR and throws away identifications (synthetic test: 2.06% estimated vs 0.62% true).

3. **PEP and q-value answer different questions; filtering at "PEP <= 0.01" is far stricter than "q <= 0.01."** PEP (posterior error probability, local FDR) is the probability that THIS PSM is wrong; q-value is the FDR of the list cut at this PSM. FDR is the average of PEP over the accepted set (Kall 2008). The worst PSM in a 1%-FDR list typically has a PEP of 10-50%. Use q-value for list cutoffs; use PEP only for per-ID decisions (e.g. picking one PTM site). And PSM-FDR at 1% does NOT give 1% peptide-FDR or 1% protein-FDR -- each level needs its own estimation; hand protein-level control to protein-inference.

## The FDR Vocabulary, Precisely

- **FDR**: the expected proportion of false positives among ALL accepted items at a threshold -- a property of the whole list.
- **q-value**: the minimum FDR at which a given PSM is still accepted; monotone after taking the running minimum from the bottom of the ranked list. Filter on q <= 0.01.
- **PEP (local FDR)**: the probability that THIS PSM is wrong given its score. Local, per-PSM; FDR is the integral of PEP over the accepted set (Kall 2008, "two sides of the same coin").
- **The estimator must match the search mode.** Concatenated target-decoy competition (TDC): FDR = (#decoy + 1) / #target (one best hit per spectrum already resolves the competition; Elias-Gygi's 2 * #decoy / (#target + #decoy) is the older, conservative whole-list form for this same composite search). Separate target and decoy searches: pi0 * #decoy / #target (Kall et al. 2008; pi0 = 1 is valid but conservative), or the refined mix-max estimator (Keich, Kertesz-Farkas & Noble 2015). Mix-max is a distinct, calibrated-score procedure for the separate-search setting.

## Tool Taxonomy

| Tool / method | Citation | Mechanism / role | When |
|---|---|---|---|
| Comet | Eng 2013 | XCorr + E-value; SEQUEST lineage, open-source | Robust default, TPP pipelines; pairs with Percolator |
| X!Tandem | -- | hyperscore + refinement passes | Legacy/free; semi-tryptic refinement niche |
| MS-GF+ | Kim & Pevzner 2014 | SpecEValue via generating-function DP | Calibrated cross-instrument E-value; ETD/CID, low-res, non-standard enzymes |
| MaxQuant / Andromeda | Cox 2011 | binomial probability score; integrated MBR/LFQ/TMT | All-in-one quant pipeline (LFQ, TMT, SILAC); GUI |
| MSFragger | Kong 2017 | hyperscore via fragment-ion indexing (~100x faster) | Open/mass-tolerant search, PTM discovery, huge datasets; core of FragPipe |
| Sage | Lazear 2023 | hyperscore-style, Rust, rescoring-native | Modern scalable open-source pipelines; emits Percolator-ready features |
| MetaMorpheus | Solntsev 2018 | calibration + G-PTM-D multinotch | PTM discovery with built-in calibration; proteoform-aware |
| pFind 3 | Chi 2018 | open-search engine | Maximal unrestricted-PTM/mutation discovery |
| Percolator | Kall 2007 | semi-supervised SVM re-rank on decoy negatives | Boost IDs at fixed FDR; non-tryptic/PTM/large search spaces |
| mokapot | Fondrie & Noble 2021 | Percolator in Python; swappable XGBoost classifier | Python pipelines, Sage output, custom features |
| MS2Rescore + DeepLC + MS2PIP | Declercq 2022; Bouwmeester 2021; Gabriels 2019 | predicted-RT + predicted-intensity rescoring features | Sharpen target/decoy separation; immunopeptidomics |
| Spectral-library search | -- | match empirical reference spectra (intensity + RT) | Faster/more specific for known peptides -> spectral-libraries |
| Protein grouping / protein FDR | Savitski 2015; The 2022 | picked / picked-group FDR | route OUT -> protein-inference |
| PTM site localization | -- | per-site PEP, localization scoring | route OUT -> ptm-analysis |

## Decision Tree by Scenario

| Scenario | Recommended | Why |
|---|---|---|
| Standard DDA, clean FDR, scriptable | Comet or Sage + Percolator/mokapot at q <= 0.01 | well-validated; rescoring boosts IDs at fixed FDR |
| Cross-instrument / varied fragmentation / odd enzyme | MS-GF+ | SpecEValue is calibrated so a threshold means the same everywhere |
| Discover unknown PTMs / mass shifts | MSFragger open search (-150..+500 Da) | fragment indexing makes wide-window search feasible; then closed search on discovered mods -> ptm-analysis |
| Huge dataset, reproducible, cloud-scale | Sage (rescoring-native) | Rust speed; emits Percolator features directly |
| All-in-one with quant in the same tool | MaxQuant/Andromeda | integrated LFQ/TMT/SILAC and MBR |
| Non-tryptic (immunopeptidomics, degradomics) | any engine + Percolator/MS2Rescore | rescoring gains are largest where search space explodes |
| Few PSMs (single-protein pulldown) | do NOT trust decoy FDR; inspect spectra manually | decoy counts too noisy below ~hundreds of PSMs |
| Need per-site / per-ID confidence | act on PEP, not q-value | q-value is list-level; PEP is local |

Default when uncertain: concatenated target-decoy search with Comet or Sage, rescore with Percolator/mokapot, filter at q <= 0.01, and hand protein-level FDR to protein-inference.

### Build the Concatenated Target-Decoy Database

**Goal:** Turn a target-only FASTA into the concatenated target+decoy database that every FDR number below depends on.

**Approach:** Reverse at the PROTEIN level with the peptide N/C termini held fixed, so decoy peptides obey the same enzyme rules and match the targets in length and amino-acid composition, at 1:1. Sage generates decoys internally and takes the TARGET-ONLY FASTA (so do MSFragger and MetaMorpheus, per their documentation); Comet, MS-GF+ and X!Tandem need this file. Include contaminants in the input before reversing, so contaminant decoys exist too. **The decoy tag is where this Skill silently fails:** every tool defaults to a different string, and a mismatch means zero decoys are found and every PSM is reported "at 1% FDR" with no error.

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
| MSFragger / FragPipe / Philosopher | `rev_` (lower case) | `decoy_prefix`, `--decoy` |
| MaxQuant | `REV__` | fixed |
| Percolator | reads the pin `Label` column | `-P` only for `--picked-protein` |

Checked on OpenMS 3.5.0: 31,437 UniProt entries (human + yeast + E. coli + contaminants) became 62,874 with 31,437 `DECOY_` entries in 5 s.

### Run a DDA Search from the Command Line

**Goal:** Search a centroided mzML against the database and emit PSMs plus Percolator features.

**Approach:** All three engines below run one concatenated search and write a Percolator `.pin`, so the rescoring step is the same for each. Convert vendor raw first (`msconvert --mzML --zlib --filter "peakPicking vendor msLevel=1-"` -> data-import). Set tolerances to the instrument, not to the default: 10 ppm precursor and high-res fragment settings for an Orbitrap, 0.6 Da / ion-trap binning for CID. `examples/dda_search.sh` runs the whole route (database -> search -> Percolator -> 1% list) for `ENGINE=sage` or `ENGINE=comet`.

```bash
# Sage 0.14.6 -- generates its own 'rev_' decoys, so it takes the TARGET-ONLY FASTA.
# search.json: {"database": {"enzyme": {"missed_cleavages": 2, "min_len": 7, "max_len": 30,
#   "cleave_at": "KR", "restrict": "P"}, "static_mods": {"C": 57.0215},
#   "variable_mods": {"M": [15.9949]}, "max_variable_mods": 2,
#   "decoy_tag": "rev_", "generate_decoys": true, "fasta": "human.fasta"},
#  "precursor_tol": {"ppm": [-10, 10]}, "fragment_tol": {"ppm": [-20, 20]},
#  "isotope_errors": [0, 1], "deisotope": true, "predict_rt": true, "report_psms": 1}
sage search.json sample.mzML --write-pin     # -> results.sage.tsv + results.sage.pin

# Comet 2026.02 -- searches the concatenated DB built above.
comet -p                                     # writes comet.params.new to edit
#   database_name = human_target_decoy.fasta   decoy_search = 0   decoy_prefix = DECOY_
#   peptide_mass_tolerance_upper = 10.0   peptide_mass_tolerance_lower = -10.0
#   peptide_mass_units = 2                 # 2 = ppm
#   fragment_bin_tol = 0.02  fragment_bin_offset = 0.0   # high-res HCD; 1.0005/0.4 for ion trap
#   search_enzyme_number = 1  allowed_missed_cleavage = 2  peptide_length_range = 7 30
#   variable_mod01 = 15.9949 M 0 2 -1 0 0 0.0   max_variable_mods_in_peptide = 2
#   output_percolatorfile = 1  output_txtfile = 1  num_output_lines = 1
comet -Pcomet.params -Nsample sample.mzML    # -> sample.pin, sample.txt, sample.pep.xml

# MS-GF+ v2024.03.26 -- calibrated SpecEValue; -tda 0 because the DB already has decoys.
java -Xmx8g -jar MSGFPlus.jar -s sample.mzML -d human_target_decoy.fasta \
  -decoy DECOY_ -o sample.mzid -t 10ppm -ti 0,1 -tda 0 \
  -m 3 -inst 3 -e 1 -ntt 2 -mod mods.txt \
  -minLength 7 -maxLength 30 -maxMissedCleavages 2 -n 1 -addFeatures 1 -thread 8
# mods.txt: "NumMods=2" / "C2H3N1O1,C,fix,any,Carbamidomethyl" / "O1,M,opt,any,Oxidation"
java -cp MSGFPlus.jar edu.ucsd.msjava.ui.MzIDToTsv -i sample.mzid -o sample.tsv -showDecoy 1
```

### Rescore to FDR-Controlled PSMs with Percolator

**Goal:** Turn engine features into a single learned score and a q-value, and read off the 1% list.

**Approach:** Percolator trains a semi-supervised SVM on the decoy PSMs with three-fold cross-validation, so it must be told which post-processing matches the search. `-Y`/`--post-processing-tdc` is target-decoy competition, correct for a concatenated search with one hit per spectrum; `-y`/`--post-processing-mix-max` is the mix-max method and **only has an effect on separate target and decoy searches**, where it is the default (flag text from `percolator --help`, 3.09.0). `--results-psms` holds targets only, so no decoy filtering afterwards -- but read the `q-value` column BY NAME, because Percolator emits a `filename` column only when the pin carries one (Sage does, Comet does not), which shifts every later column by one.

```bash
percolator --post-processing-tdc \
  --results-psms psms.target.tsv       --decoy-results-psms psms.decoy.tsv \
  --results-peptides peptides.target.tsv --decoy-results-peptides peptides.decoy.tsv \
  results.sage.pin

# 1% list, with the q-value column located by name rather than by index
awk -F'\t' 'NR == 1 { for (i = 1; i <= NF; i++) if ($i == "q-value") q = i; next }
            q && $q <= 0.01' psms.target.tsv > psms_1pct.tsv
```

Rescoring pays where the engine's own score is weakest. On one Orbitrap Astral 5-min DDA run (250 pg HYE load, 6,135 MS2 spectra, 31,437-protein database, PXD070049), Comet's raw `-log10(e-value)` gave **916** PSMs at 1% FDR and Percolator lifted the same search to **1,144** (+25%); Sage's `sage_discriminant_score` is already a learned score, so Percolator moved it from **1,406** to **1,398** (-0.6%) -- no gain to be had. Sage 0.14.6 wins this input outright; Comet 2026.02 with Percolator lands 19% behind it, and MS-GF+ v2024.03.26 read straight off `-log10(SpecEValue)` with no rescoring gives **656**, because a calibrated E-value buys cross-instrument comparability, not raw yield. **Do not generalise these counts**: one run, one low-load short-gradient method, each engine at its own idiomatic high-res settings. Rescoring also needs training data -- pooling the three DDA runs (6,864 PSMs) gave Percolator 5,006 PSMs against Sage's own 4,961, while on a single run of 1,939 PSMs it had too few positives to improve anything.

### Database Search with pyOpenMS

**Goal:** Match tandem mass spectra in an mzML file against a protein FASTA and produce scored PSMs as idXML.

**Approach:** `SimpleSearchEngineAlgorithm` actually scores spectra (the hand-rolled `ProteaseDigestion` loop only digests, it never matches a spectrum). The FASTA must already contain target + decoy sequences concatenated for downstream FDR; decoys carry a recognizable prefix (or set `decoys` to `'true'` to let the engine generate them). Defaults are 10 ppm fragment tolerance and 1 missed cleavage, so set parameters explicitly. The search already annotates `target_decoy` on each hit.

```python
from pyopenms import SimpleSearchEngineAlgorithm, IdXMLFile, PeptideIdentificationList

protein_ids = []
peptide_ids = PeptideIdentificationList()   # pyOpenMS 3.5+: a plain [] raises TypeError
search = SimpleSearchEngineAlgorithm()
p = search.getParameters()
p.setValue('precursor:mass_tolerance', 10.0)
p.setValue('precursor:mass_tolerance_unit', 'ppm')
p.setValue('fragment:mass_tolerance', 0.02)         # HCD Orbitrap; default is 10 ppm
p.setValue('fragment:mass_tolerance_unit', 'Da')
p.setValue('peptide:missed_cleavages', 2)           # default is 1
p.setValue('modifications:fixed', [b'Carbamidomethyl (C)'])
p.setValue('modifications:variable', [b'Oxidation (M)'])
search.setParameters(p)
# spectra are scored against in-silico fragment ions of every candidate peptide
search.search('sample.mzML', 'human_target_decoy.fasta', protein_ids, peptide_ids)

# protein_ids FIRST in load/store -- the OpenMS argument order is fixed
IdXMLFile().store('search_results.idXML', protein_ids, peptide_ids)
```

### Annotate Target/Decoy and Estimate FDR with pyOpenMS

**Goal:** Convert raw PSM scores into q-values and keep only PSMs at 1% FDR.

**Approach:** `PeptideIndexing` maps each PSM back to proteins and flags target vs decoy from the decoy prefix -- needed for idXML from other engines or after changing the FASTA; `SimpleSearchEngineAlgorithm` output above is already annotated and can go straight to `FalseDiscoveryRate`. `FalseDiscoveryRate.apply` runs the concatenated competition; `IDFilter` keeps q <= 0.01. This is the real pyOpenMS path -- not a hand-rolled decoy/target ratio of unknown provenance.

```python
from pyopenms import PeptideIndexing, FalseDiscoveryRate, IDFilter, FASTAFile

fasta = []
FASTAFile().load('human_target_decoy.fasta', fasta)
indexer = PeptideIndexing()
params = indexer.getParameters()
params.setValue('decoy_string', 'DECOY_')      # must match the decoy prefix in the FASTA
params.setValue('decoy_string_position', 'prefix')
indexer.setParameters(params)
indexer.run(fasta, protein_ids, peptide_ids)   # sets target/decoy flags on every hit

FalseDiscoveryRate().apply(peptide_ids)         # concatenated competition -> per-PSM q-value as the new score
IDFilter().filterHitsByScore(peptide_ids, 0.01) # 0.01 = 1% FDR, the community list-level standard
IDFilter().removeDecoyHits(peptide_ids)
```

### FDR from a Results Table (concatenated competition, made explicit)

**Goal:** Compute q-values from any engine's PSM table when the search was a single concatenated target-decoy search.

**Approach:** Keep the best hit per spectrum, rank by score, walk down accumulating target and decoy counts, FDR = (decoys + 1)/targets, then take the running minimum from the bottom to get monotone q-values. `score` must be higher-is-better, and the decoy prefix must match the engine's (Sage and FragPipe write lowercase `rev_`); a table with no recognised decoys must stop, not pass every PSM. This form is correct ONLY for concatenated competition; separate searches need pi0 * decoys/targets (Kall et al. 2008; pi0 = 1 is the conservative default) or the mix-max estimator (Keich, Kertesz-Farkas & Noble 2015; Percolator's default for separate-search input).

```python
import pandas as pd

DECOY_PREFIXES = ('decoy_', 'rev_', 'xxx_')   # compared lower-cased: DECOY_, REV_ / REV__ (MaxQuant), rev_ (Sage, FragPipe), XXX_

psms = pd.read_csv('search_results.tsv', sep='\t')   # map engine columns to 'scan', 'score', 'protein'
# score must be HIGHER-is-better: Sage sage_discriminant_score or Comet xcorr as is;
# E-values (Comet e-value, MS-GF+ SpecEValue) as -log10(E-value)
psms['is_decoy'] = psms['protein'].str.lower().str.startswith(DECOY_PREFIXES)
if not psms['is_decoy'].any():
    raise ValueError('no decoy PSMs recognised: check the decoy prefix, or the table was already decoy-filtered')
# one best hit per spectrum (Comet .txt writes 5 rows per scan by default)
psms = psms.sort_values('score', ascending=False).drop_duplicates('scan').reset_index(drop=True)

# concatenated target-decoy competition: each decoy above threshold estimates one false target
targets = (~psms['is_decoy']).cumsum()
decoys = psms['is_decoy'].cumsum()
psms['fdr'] = (decoys + 1) / targets.clip(lower=1)   # +1: zero decoys is not zero FDR (OpenMS conservative default)
psms['qvalue'] = psms['fdr'][::-1].cummin()[::-1]   # running min from the bottom -> monotone q-values

kept = psms[(psms['qvalue'] <= 0.01) & (~psms['is_decoy'])]   # 1% list-level FDR
```

### FDR from SEPARATE Target and Decoy Searches (pi0 * D / T)

**Goal:** Get a valid 1% list out of two result tables produced by searching the same spectra against a target DB and a decoy DB independently.

**Approach:** No competition resolved which hit wins, so the decoy count estimates the number of incorrect TARGETS directly, scaled by pi0, the proportion of target PSMs that are incorrect (Kall et al. 2008). pi0 = 1 is always valid and conservative; the median-decoy estimate (twice the fraction of target scores below the median decoy score) recovers the identifications that pi0 = 1 throws away, at the cost of estimating a nuisance parameter. Do NOT run the concatenated snippet above on the two tables merged. `examples/separate_search_fdr.py` ships this as a script.

```python
import numpy as np

def estimate_pi0(target_scores, decoy_scores):        # Kall et al. 2008
    median_decoy = np.median(decoy_scores)
    return min(1.0, 2.0 * np.mean(np.asarray(target_scores) < median_decoy))

def separate_search_qvalues(targets, decoys, pi0=None):
    # targets, decoys: DataFrames with 'scan' and 'score', ONE row per spectrum each
    if pi0 is None:
        pi0 = estimate_pi0(targets['score'].to_numpy(), decoys['score'].to_numpy())
    t = targets.sort_values('score', ascending=False).reset_index(drop=True)
    decoy_sorted = np.sort(decoys['score'].to_numpy())
    n_decoy_above = len(decoy_sorted) - np.searchsorted(decoy_sorted, t['score'].to_numpy(), side='left')
    t['fdr'] = np.minimum(1.0, pi0 * n_decoy_above / np.arange(1, len(t) + 1))
    t['qvalue'] = t['fdr'][::-1].cummin()[::-1]
    return t, pi0
```

On the synthetic separate-search pair with ground truth (12,000 spectra each): pi0-hat = 0.611 keeps 2,888 PSMs at a true FDP of 1.04%, pi0 = 1 keeps 2,588 at 0.62%, and Elias-Gygi's 2d/(t+d) misapplied here keeps only 2,139 at 0.33%. The alternative is to hand both tables to Percolator and let mix-max do it: that is Percolator's default for separate-search input, and it is a calibrated-score procedure, not the same arithmetic.

## Per-Method Failure Modes

### Concatenated vs separate FDR formula mismatch
**Trigger:** running the concatenated snippet on a merged table from separate searches without per-spectrum competition, or applying Elias-Gygi's 2*decoy/(target+decoy) to separate searches.
**Mechanism:** Elias-Gygi's factor 2 counts decoys in the combined target+decoy list of a concatenated search; in separate searches every spectrum gets both a target and a decoy hit, so the decoy count estimates false targets directly, scaled by pi0 (the fraction of target PSMs that are incorrect).
**Symptom:** mis-estimated FDR; 2d/(t+d) on separate searches over-estimates it (synthetic test: 2.06% vs 0.62% true, about a quarter of IDs lost).
**Fix:** confirm the search mode; concatenated -> (#decoy + 1)/#target; separate -> pi0 * #decoy/#target (Kall et al. 2008; code in "FDR from SEPARATE Target and Decoy Searches" above and in `examples/separate_search_fdr.py`) or the mix-max estimator (Keich, Kertesz-Farkas & Noble 2015). In Percolator, mix-max is the default for separate-search input and `-Y`/`--post-processing-tdc` selects target-decoy competition instead; concatenated input forces TDC automatically.

### Thresholding on raw engine score
**Trigger:** filtering on XCorr/hyperscore/Andromeda score, or comparing scores from two engines.
**Mechanism:** scores are uncalibrated, charge/length-dependent, and not monotone in true probability.
**Symptom:** different cutoffs admit different real FDRs; cross-engine merges nonsensical.
**Fix:** always convert to q-value (or SpecEValue/PEP) first; rescore with Percolator/mokapot.

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
**Fix:** use Percolator/mokapot default cross-validation; predicted-feature rescoring (DeepLC/MS2PIP) is safer; validate with entrapment for high-stakes claims (Wen 2025).

### DIA tool FDR taken at face value
**Trigger:** trusting a DIA tool's reported 1% peptide/protein FDR.
**Mechanism:** entrapment shows several DIA tools do not reliably control FDR (Wen 2025).
**Symptom:** real error rate exceeds the reported FDR.
**Fix:** validate with entrapment for high-stakes DIA claims -> dia-analysis.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|---|---|---|
| Precursor tolerance 10-20 ppm (high-res Orbitrap) | -- | matches FT mass accuracy; tighter = fewer random candidates at fixed FDR |
| Precursor tolerance -150..+500 Da (open search) | Kong 2017 | captures arbitrary PTM/mutation shifts; feasible only with fragment indexing |
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
| Every PSM passes 1% FDR, or `ValueError: no decoy PSMs recognised` | decoy prefix not matched (Sage and FragPipe write lowercase `rev_`) | compare prefixes lower-cased; check `psms['protein'].str[:6].value_counts()` |
| Empty 1% list from a table snippet | lower-is-better score (E-value, SpecEValue) used as `score` | use `-log10(E-value)` |
| All q-values 0 from a hand-rolled table | no +1 correction on a list with zero decoys | use (decoys + 1)/targets; a tiny list cannot reach 1% |
| Percolator q-method mismatched to search mode | mix-max is the default for separate-search input | for separate searches, mix-max (default) or `-Y`/`--post-processing-tdc` for target-decoy competition; concatenated input forces TDC automatically; use `--picked-protein` for protein FDR |
| Percolator's 1% list is far too big or too small when cut by column index | Percolator writes a `filename` column only when the pin has one (Sage yes, Comet no), shifting `q-value` between columns 3 and 4 | locate `q-value` by header name, never by a fixed index |
| Decoy count is twice the target count after a search | `-tda 1` (MS-GF+) or `generate_decoys: true` (Sage) run against a database that already contains decoys | use `-tda 0` / `generate_decoys: false` with a concatenated DB, or feed the target-only FASTA and let the engine make them |
| `philosopher peptideprophet` exits 0 but the output has no `peptideprophet_result` and the log says `read in 0 1+, 0 2+ ... spectra` / `read in no data` | its embedded PeptideProphet does not model this pepXML (seen with Comet 2026.02 rev.2 pepXML under Philosopher v5.1.0 on Windows, with and without `--nonparam --decoy`) | rescore with Percolator on the engine's `.pin` instead; check the interact file for `peptideprophet_result` before trusting a PeptideProphet run |
| 1% PSM FDR assumed to give 1% protein FDR | each level needs its own estimation | estimate protein-level (picked) FDR -> protein-inference |
| "PEP <= 0.01" returns far fewer IDs than expected | PEP is per-PSM and far stricter than q-value | filter list cutoffs on q-value; reserve PEP for per-ID decisions |

## References

- Elias, J.E. & Gygi, S.P. 2007. Target-decoy search strategy for increased confidence in large-scale protein identifications by mass spectrometry. *Nature Methods* 4(3):207-214.
- Keich, U., Kertesz-Farkas, A. & Noble, W.S. 2015. Improved false discovery rate estimation procedure for shotgun proteomics. *Journal of Proteome Research* 14(8):3148-3161.
- Kall, L., Canterbury, J.D., Weston, J., Noble, W.S. & MacCoss, M.J. 2007. Semi-supervised learning for peptide identification from shotgun proteomics datasets. *Nature Methods* 4(11):923-925.
- Kall, L., Storey, J.D., MacCoss, M.J. & Noble, W.S. 2008. Assigning significance to peptides identified by tandem mass spectrometry using decoy databases. *Journal of Proteome Research* 7(1):29-34.
- Kall, L., Storey, J.D., MacCoss, M.J. & Noble, W.S. 2008. Posterior error probabilities and false discovery rates: two sides of the same coin. *Journal of Proteome Research* 7(1):40-44.
- Eng, J.K., Jahan, T.A. & Hoopmann, M.R. 2013. Comet: an open-source MS/MS sequence database search tool. *Proteomics* 13(1):22-24.
- Kim, S. & Pevzner, P.A. 2014. MS-GF+ makes progress towards a universal database search tool for proteomics. *Nature Communications* 5:5277.
- Cox, J., Neuhauser, N., Michalski, A., Scheltema, R.A., Olsen, J.V. & Mann, M. 2011. Andromeda: a peptide search engine integrated into the MaxQuant environment. *Journal of Proteome Research* 10(4):1794-1805.
- Kong, A.T., Leprevost, F.V., Avtonomov, D.M., Mellacheruvu, D. & Nesvizhskii, A.I. 2017. MSFragger: ultrafast and comprehensive peptide identification in mass spectrometry-based proteomics. *Nature Methods* 14(5):513-520.
- Lazear, M.R. 2023. Sage: an open-source tool for fast proteomics searching and quantification at scale. *Journal of Proteome Research* 22(11):3652-3659.
- Solntsev, S.K., Shortreed, M.R., Frey, B.L. & Smith, L.M. 2018. Enhanced global post-translational modification discovery with MetaMorpheus. *Journal of Proteome Research* 17(5):1844-1851.
- Chi, H., Liu, C., Yang, H. et al. 2018. Comprehensive identification of peptides in tandem mass spectra using an efficient open search engine. *Nature Biotechnology* 36:1059-1061.
- Fondrie, W.E. & Noble, W.S. 2021. mokapot: fast and flexible semisupervised learning for peptide detection. *Journal of Proteome Research* 20(4):1966-1971.
- Bouwmeester, R., Gabriels, R., Hulstaert, N., Martens, L. & Degroeve, S. 2021. DeepLC can predict retention times for peptides that carry as-yet unseen modifications. *Nature Methods* 18:1363-1369.
- Gabriels, R., Martens, L. & Degroeve, S. 2019. Updated MS2PIP web server delivers fast and accurate MS2 peak intensity prediction for multiple fragmentation methods, instruments and labeling techniques. *Nucleic Acids Research* 47(W1):W295-W299.
- Declercq, A., Bouwmeester, R., Hirschler, A., Carapito, C., Degroeve, S., Martens, L. & Gabriels, R. 2022. MS2Rescore: data-driven rescoring dramatically boosts immunopeptide identification rates. *Molecular & Cellular Proteomics* 21(8):100266.
- Savitski, M.M., Wilhelm, M., Hahne, H., Kuster, B. & Bantscheff, M. 2015. A scalable approach for protein false discovery rate estimation in large proteomic data sets. *Molecular & Cellular Proteomics* 14(9):2394-2404.
- The, M., Samaras, P., Kuster, B. & Wilhelm, M. 2022. Reanalysis of ProteomicsDB using an accurate, sensitive, and scalable false discovery rate estimation approach for protein groups. *Molecular & Cellular Proteomics* 21(12):100437.
- Wen, B., Freestone, J., Riffle, M., MacCoss, M.J., Noble, W.S. & Keich, U. 2025. Assessment of false discovery rate control in tandem mass spectrometry analysis using entrapment. *Nature Methods* 22:1454-1463.

## Related Skills

- protein-inference - Group peptides to protein groups and control protein-level (picked) FDR
- ptm-analysis - Open/variable-mod search follow-up and per-site PTM localization
- dia-analysis - DIA peptide-centric extraction and scoring; entrapment FDR validation
- quantification - FDR-filtered IDs feed label-free/TMT intensity quantification
- spectral-libraries - Empirical and predicted spectral-library search as an ID alternative
- data-import - Load mzML/raw MS data before identification
- database-access/uniprot-access - Build the target FASTA (canonical vs isoform, contaminants)
