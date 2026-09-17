---
name: bio-proteomics-dia-analysis
description: Analyzes data-independent acquisition (DIA) proteomics by scoring reconstructed fragment-chromatogram peak groups against a decoy null with DIA-NN (library-free search against a deep-learning predicted library, or library-based search), Spectronaut, OpenSWATH, and EncyclopeDIA. Frames the deliverable around q-value LEVEL (precursor/peptide/protein-group) and CONTEXT (run vs experiment-wide/global) rather than a bare "1% FDR", and around the duty-cycle-vs-selectivity acquisition tradeoff (window design, staggered demultiplexing, diaPASEF, narrow-window Astral). Use when identifying and quantifying proteins from DIA mass spectrometry runs and filtering DIA-NN report.parquet/matrix output. Building the spectral library itself is spectral-libraries; normalization and protein roll-up is quantification; statistical testing of the matrix is differential-abundance.
tool_type: cli
primary_tool: diann
license: MIT
---

## Version Compatibility

Reference examples tested with: pandas 3.0.5 (report filter, run); DIA-NN flags checked against the DIA-NN README at 1.9.2, 2.0 and master 2.6.1 (not executed); easypqp 0.1.59 CLI source (checked 2026-09-15)

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# DIA Analysis -- Scoring Reconstructed Peak Groups Against a Decoy Null, Filtered at the Right q-Value Context

**"Identify and quantify proteins from my DIA runs"** -> Reconstruct, per candidate peptide, a set of co-eluting fragment extracted-ion chromatograms (XICs) and score whether that peak group is real against a decoy null -- because every wide-isolation-window MS2 is chimeric, so the problem is deconvolution and peak-group scoring, not spectrum matching.
- CLI: `diann --fasta-search --predictor` for library-free search (in DIA-NN this IS the predicted-library route; the modern default)
- CLI: `diann --lib predicted.speclib` to reuse an in-silico library predicted earlier
- CLI: `diann --lib empirical.parquet` for an experimental or chromatogram library
- CLI: `OpenSwathWorkflow` + `pyprophet` when explicit run/experiment/global FDR contexts must be auditable

Scope: this skill OWNS running the DIA search engine and FILTERING its output at the correct q-value level and context. Building the library (experimental, chromatogram, predicted) -> spectral-libraries. Normalization, MaxLFQ roll-up, and matrix summarization -> quantification. Statistical testing of the protein matrix -> differential-abundance. Loading raw vendor/mzML data -> data-import. OUT OF SCOPE: DDA spectrum-to-peptide matching (peptide-identification); pathway enrichment of the hit list; acquiring the data (the analyst inherits the window design from the core facility).

## The Single Most Important Modern Insight -- The q-Value Context Is a Study-Design Choice, Not a Default

1. **DIA quantification is peak-group SCORING against a decoy null, not spectrum matching.** Gillet 2012 inverted DDA: instead of asking "what peptide is this spectrum", DIA asks, per library peptide, "does a co-eluting peak group of this peptide's expected fragments exist in the chimeric MS2 stream". Decoys are shuffled or reversed peptide queries scored identically; the q-value is the expected fraction of accepted IDs that are decoy-like false peak groups. The count of "proteins found" is therefore a function of the decoy-calibrated threshold, never a quality metric in itself.

2. **"1% FDR" is meaningless without naming the LEVEL and the CONTEXT -- state both.** LEVEL = precursor vs peptide vs protein-group (filtering precursors at 1% does NOT give proteins at 1%; control both). CONTEXT = run-specific vs experiment-wide vs global (Rosenberger 2017). Naively filtering N runs at per-run 1% inflates the experiment-wide error: 1% per run accumulates false positives across the union, severe at hundreds-to-thousands of runs. For a cross-run matrix, filter on the GLOBAL protein-group q-value, not the per-run one. The column chosen (`Q.Value` vs `Global.PG.Q.Value`) is the decision.

3. **The predicted-library route is now the default recommendation.** Library-based search is sensitive but capped by an ill-matched library (wrong organism/tissue/mods silently limits coverage with no error). Library-free directDIA finds sample-specific content but its larger implicit search space can INFLATE IDs if FDR is not controlled across the two-pass process -- worst on wide-window chimeric data. The compromise the field converged on: predict an in-silico library for the whole FASTA digest (DIA-NN built-in predictor, or Prosit/AlphaPeptDeep) and search against THAT, getting directDIA's "no wet-lab library" with library-based's bounded, better-calibrated search. In DIA-NN these are not separate routes: its "library-free" search is exactly this search against a library predicted from the whole FASTA digest.

## Chimerism, Deconvolution, and the Window Design the Analyst Inherits

DDA picks top-N precursors and fragments each in isolation, so every MS2 is nominally one peptide. DIA abandons selection: the quadrupole steps through wide isolation windows (4-25 Th classically, 2 Th on Astral, mobility-gated on timsTOF) and co-fragments EVERY precursor in each window. Consequence: every DIA MS2 is chimeric, a superposition of fragments from all co-isolated precursors. The engine must deconvolve -- reconstruct each candidate's fragment XICs and score the peak group.

Selectivity is set by isolation-window WIDTH. Narrower window = fewer co-isolated precursors = less chimerism = cleaner XICs = fewer false peak groups. But narrower windows mean MORE windows per cycle -> longer duty cycle -> fewer points across each LC peak -> worse quant precision. The central acquisition tradeoff is duty cycle (sampling speed) vs selectivity (window width). Rule of thumb: aim for >= 6 MS2 points across the FWHM of an LC peak for reliable quant. The analyst INHERITS this design and must not pretend all DIA is equivalent:

- Fixed windows (classic SWATH): 32 x 25 Th. Simple; wastes selectivity because precursor density is non-uniform across m/z.
- Variable windows: widths chosen so each holds roughly equal precursor density (narrow where the proteome is dense ~600-800 m/z). Orbitrap best practice.
- Staggered / overlapping windows + demultiplexing (Amodei 2019): two interleaved patterns offset by half a window; demultiplexing recovers effective windows of half the physical width without halving duty cycle. DIA-NN supports overlapping-window schemes natively (DIA-NN docs), and skipping demultiplexing cost only a few % coverage in a community benchmark (DIA-NN GitHub #438), so demultiplexing is optional there. Engines without native support need demultiplexed input, converted with vendor peak picking as the FIRST filter (msconvert warns that vendor peak picking has no effect otherwise): `msconvert run.raw --mzML --filter "peakPicking vendor msLevel=1-" --filter "demultiplex optimization=overlap_only massError=10ppm"`. Current DIA-NN also reads Thermo .raw directly, including native Linux builds. MSX (randomized window combinations) is largely historical.
- diaPASEF (Meier 2020): on timsTOF the isolation tile rides the m/z-vs-ion-mobility diagonal, so only precursors sharing BOTH m/z AND mobility co-isolate -- the mobility dimension is an orthogonal selectivity filter for free. The engine extracts a 4D peak group (RT x m/z x fragment x mobility).
- Narrow-window Astral (Guzman 2024): >200 Hz MS/MS makes 2-Th windows feasible across the whole range; at 2 Th the MS2 is nearly non-chimeric, which makes library-free directDIA far more trustworthy than it was on 25-Th SWATH.

## Tool Taxonomy

| Tool / method | Citation | Mechanism / role | When |
|---------------|----------|------------------|------|
| DIA-NN | Demichev 2020 | Deep-NN peak-group scoring + interference correction + QuantUMS quant; library-free, predicted, or library-based | Default for high-throughput, large cohorts, diaPASEF, Astral; free, scriptable |
| Spectronaut | Biognosys (commercial) | directDIA+ pipeline with in-app DL prediction; mature GUI/QC | Regulated/clinical work, polished QC, mixed vendors, when licensed |
| OpenSWATH + PyProphet | Rost 2014; Rosenberger 2017 | Classic peptide-centric extraction + semi-supervised scoring with explicit run/experiment/global q-contexts | When auditable FDR-context control is required; library-based ONLY, needs iRT/RT alignment |
| EncyclopeDIA / Walnut | Searle 2018 | Chromatogram-library search (.dlib/.elib) + GPF; Walnut = library-free mode | Building project-specific chromatogram-library depth on Orbitrap |
| FragPipe (MSFragger-DIA / DIA-Umpire) | -- | Spectrum-centric via pseudo-spectra + Philosopher FDR; IonQuant explicit MBR-FDR | Unified DDA+DIA shop in the MSFragger ecosystem |
| Skyline | MacCoss lab | Targeted/visual peak inspection and demultiplexing; not a discovery engine | Manual peak curation, PRM, library curation -- the microscope, not the engine |
| AlphaDIA | Mann lab | Open transformer-based end-to-end Python; AlphaPeptDeep predictions | Cutting-edge open research, Astral, Python-native pipelines |
| Library build | (route OUT) | Experimental/chromatogram/predicted library construction | -> spectral-libraries |
| Stats on the matrix | (route OUT) | Normalization, roll-up, moderated testing | -> quantification, differential-abundance |

Tool leadership moves fast (Astral, AlphaDIA, DIA-NN releases). Confirm the current recommended engine and version for the specific instrument before committing rather than hard-coding "DIA-NN is best".

## Decision Tree by Scenario

| Scenario | Recommended | Why |
|----------|-------------|-----|
| Discovery cohort, no wet-lab library | `diann --fasta-search --gen-spec-lib` (predicted route) then search against it | Bounded, better-calibrated search vs raw directDIA; the modern default |
| Quick single-run discovery, Astral 2-Th data | DIA-NN library-free (directDIA) | Near-non-chimeric MS2 makes directDIA trustworthy |
| Have a deep experimental/chromatogram library | `diann --lib library.parquet --fasta db.fasta` (no `--fasta-search`) | Targeted extraction is most sensitive when the library matches |
| Need auditable run/experiment/global FDR for a regulated submission | OpenSWATH + PyProphet | Explicit q-value contexts per Rosenberger 2017 |
| Building chromatogram-library depth for one project on Orbitrap | EncyclopeDIA (GPF) -> .elib -> DIA-NN | Empirical RT and real fragmentation in the project's own LC |
| Large cohort (hundreds-thousands of runs) | Empirical library from 20-100 high-quality runs (library-free search + `--gen-spec-lib`), then search ALL runs with that library, MBR off, fixed mass accuracies; filter on `Global.Q.Value`/`Global.PG.Q.Value` | DIA-NN docs: MBR is a convenience feature (second pass held in memory); one fixed library and fixed tolerances keep batches comparable |
| Staggered/overlapping acquisition | DIA-NN directly (native support); for engines without it, demultiplex at conversion with peak picking first | Demultiplexing is an optional few-% gain in DIA-NN; required where the engine assumes non-overlapping windows |
| PTM / peptidoform-resolved work | DIA-NN `--peptidoforms` + matched variable mods | Peptidoform-resolved target-decoy scoring |

Default when uncertain: DIA-NN with the predicted-library route (`--fasta-search --gen-spec-lib --reanalyse`), mass accuracies fixed for the instrument, then filter `Q.Value`, `PG.Q.Value`, `Global.Q.Value` and `Global.PG.Q.Value` at 0.01 for cross-run matrices (DIA-NN 1.9.x with MBR: `Lib.Q.Value`/`Lib.PG.Q.Value` instead of the `Global.*` pair).

## DIA-NN -- Predicted-Library (directDIA) Route

The default route: digest the FASTA in silico, predict a library, and search the DIA data against it in one command. Mass accuracies of 0 (the default) mean automatic: DIA-NN optimises them on the FIRST run and reuses them for every other run, so results depend on run order (`--individual-mass-acc` optimises per run). For production analyses the DIA-NN docs prefer fixed values per instrument: timsTOF `--mass-acc 15 --mass-acc-ms1 15` (shown below), Orbitrap Astral `--mass-acc 10 --mass-acc-ms1 4`, TripleTOF 6600/ZenoTOF `--mass-acc 20 --mass-acc-ms1 20`. `--reanalyse` enables MBR: a second pass using the empirical library built in the first. `--qvalue 0.01` tightens DIA-NN's default 5% precursor filter.

```bash
diann \
    --f sample1.mzML --f sample2.mzML \
    --lib "" --fasta uniprot_human.fasta --fasta-search \
    --gen-spec-lib --predictor \
    --out diann_out/report.parquet \
    --out-lib diann_out/report-lib.parquet \
    --qvalue 0.01 \
    --matrices \
    --mass-acc 15 --mass-acc-ms1 15 \
    --reanalyse --smart-profiling \
    --cut K*,R* --missed-cleavages 1 \
    --min-pep-len 7 --max-pep-len 30 \
    --unimod4 --var-mods 1 --var-mod UniMod:35,15.994915,M \
    --threads 8
```

## DIA-NN -- Library-Based Route

Supply an existing library (experimental, chromatogram-derived, or a previously predicted `.speclib`). Omit `--fasta-search`: extraction is targeted to the library content. Still pass `--fasta` so proteins are annotated from the database (add `--reannotate` for a third-party library). DIA-NN recommends DIA-NN-generated `.parquet`/`.speclib` libraries; for third-party libraries (e.g. FragPipe) compatibility should be verified for each software version, settings and data type.

```bash
diann \
    --f sample1.mzML --f sample2.mzML \
    --lib spectral_library.parquet --fasta uniprot_human.fasta \
    --out diann_out/report.parquet \
    --qvalue 0.01 --matrices \
    --mass-acc 15 --mass-acc-ms1 15 \
    --reanalyse --smart-profiling \
    --threads 8
```

## DIA-NN Output and Correct Filtering

DIA-NN 1.9+ writes the main report as Apache Parquet (`report.parquet`) by default; 2.0 makes it the only default. Matrices stay TSV. Pipelines hard-coding `report.tsv` silently break or read a stale file -- read parquet. Filter on q-value columns BEFORE pivoting to a matrix.

```
report.parquet          # main report (1.9+ default; was report.tsv pre-1.9)
report.stats.tsv        # per-run statistics
report.pg_matrix.tsv    # protein-group wide matrix
report.pr_matrix.tsv    # precursor wide matrix (verify exact dotting vs installed version)
report.gg_matrix.tsv    # gene-group wide matrix
report-lib.parquet      # generated empirical library (.parquet since 1.9.1; predicted libraries are *.predicted.speclib)
```

```python
import pandas as pd, numpy as np
report = pd.read_parquet('diann_out/report.parquet')  # NOT report.tsv on 1.9+

# Run-level LEVELS (Q.Value, PG.Q.Value) plus experiment-wide CONTEXT (Global.Q.Value, Global.PG.Q.Value) for a
# cross-run matrix. DIA-NN 1.9.x with MBR: use Lib.Q.Value / Lib.PG.Q.Value instead of the Global.* pair.
filt = report[(report['Q.Value'] <= 0.01) &
              (report['PG.Q.Value'] <= 0.01) &
              (report['Global.Q.Value'] <= 0.01) &
              (report['Global.PG.Q.Value'] <= 0.01)]  # 0.01 = standard 1% FDR

# Pivot to a protein matrix from the filtered long report.
pg = filt.pivot_table(index='Protein.Group', columns='Run', values='PG.MaxLFQ', aggfunc='first')
pg = np.log2(pg.replace(0, np.nan))  # DIA-NN writes 0 for not-quantified; log2(0) = -Inf
```

The `*_matrix.tsv` files apply an EXTRA 5% run-specific protein FDR (DIA-NN default, `--matrix-spec-q`), so the matrix protein count can be lower than the report count. A report-vs-matrix mismatch is EXPECTED, not a bug -- do not panic and do not compare the two counts as if they should match.

## Building an Empirical Library from FragPipe DDA Results (EasyPQP)

DIA-NN's built-in predictor covers predicted libraries. To build an EMPIRICAL library from FragPipe DDA search results, EasyPQP first converts each run's pepXML and spectra into pickles, then assembles the library from all of them; FragPipe can also emit a DIA-NN-format library directly.

```bash
# Real easypqp subcommands: convert, library, insilico-library (NOT a convert --format diann).
# 1) per run: pepXML + spectra (mzXML; MGF for timsTOF) -> <run>.psmpkl and <run>.peakpkl ('psmpkl'/'peakpkl' in the names)
easypqp convert --pepxml interact-run1.pep.xml --spectra run1.mzXML --psms run1.psmpkl --peaks run1.peakpkl
# 2) library from all pickles; --psmtsv requires --peptidetsv
easypqp library \
    --psmtsv psm.tsv --peptidetsv peptide.tsv \
    --rt_reference irt.tsv \
    --out library.tsv \
    *.psmpkl *.peakpkl
# With psm.tsv + peptide.tsv given, easypqp ignores --psm/--peptide/--protein_fdr_threshold (FragPipe already filtered).
# FragPipe's DIA workflow can emit a DIA-NN-format library directly -- prefer that when in FragPipe.
```

## Per-Method Failure Modes

### Library-free (directDIA)
**Trigger:** `--fasta-search` on wide-window (25-Th SWATH) data without two-pass global FDR.
**Mechanism:** building the library from the same data then quantifying it reuses the data twice; the large implicit search space inflates IDs when decoys are not controlled across both passes.
**Symptom:** implausibly high protein counts, poor reproducibility across replicates.
**Fix:** keep `--reanalyse` (two-pass global FDR) ON and filter on `Global.PG.Q.Value`; prefer the predicted-library route on chimeric data.

### Library-based
**Trigger:** library organism/tissue/modifications/gradient do not match the sample.
**Mechanism:** targeted extraction can only find what is in the library; a mismatched library caps coverage with no error.
**Symptom:** low ID counts, conserved-peptide bias (e.g. a human library on a mouse sample finds only conserved peptides).
**Fix:** match the library to the biology (mods especially); regenerate via spectral-libraries or switch to the predicted route.

### Cross-run cohort FDR
**Trigger:** filtering N runs at per-run `Q.Value <= 0.01` and unioning.
**Mechanism:** 1% per run accumulates across the union -> experiment-wide error far above 1%.
**Symptom:** inflated total protein list; irreproducible "hits" in differential testing.
**Fix:** filter `Global.Q.Value <= 0.01` and `Global.PG.Q.Value <= 0.01` for the matrix; on DIA-NN 1.9.x with MBR use `Lib.Q.Value`/`Lib.PG.Q.Value` instead (1.9.2 docs; the note is absent from 2.x docs). For very large cohorts, search all runs with an empirical library built from 20-100 runs.

### Match-between-runs (MBR)
**Trigger:** treating MBR-transferred quant values as equally confident as directly identified ones.
**Mechanism:** transferring an ID by RT/m/z/CCS matching can be a false transfer, especially for low-abundance precursors; MBR has its OWN FDR.
**Symptom:** spurious low-abundance quant filling missing values that should stay missing.
**Fix:** rely on DIA-NN's global/empirical-library q-values (IonQuant uses an explicit MBR-FDR model); do not disable `--reanalyse` then trust per-run counts.

### Staggered acquisition
**Trigger:** running staggered/overlapping data without demultiplexing through an engine that lacks native overlapping-window support.
**Mechanism:** the engine sees the wide physical window, keeping all the interference the staggering was meant to remove.
**Symptom:** noisy directDIA, poor selectivity on data that should be clean.
**Fix:** demultiplex at conversion with peak picking first (`--filter "peakPicking vendor msLevel=1-" --filter "demultiplex optimization=overlap_only massError=10ppm"`), or use DIA-NN, which handles overlapping windows natively.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| Precursor q-value `<= 0.01` | `--qvalue 0.01` (DIA-NN default is 0.05) | Standard 1% per-precursor FDR (run context); the main report is otherwise filtered only at 5%. |
| `Global.Q.Value <= 0.01` | DIA-NN docs | Experiment-wide precursor FDR, recommended alongside `Global.PG.Q.Value`. |
| `Global.PG.Q.Value <= 0.01` | Rosenberger 2017 | Experiment-wide protein-group FDR for cross-run matrices; run-specific PG q-value is NOT enough for a cohort. |
| Matrix run-specific PG filter `0.05` | DIA-NN default (`--matrix-spec-q`) | Extra 5% run-specific protein FDR applied only when building matrices -> report vs matrix count differs (expected). |
| `Lib.Q.Value`, `Lib.PG.Q.Value <= 0.01` | DIA-NN 1.9.2 docs | With MBR on 1.9.x these replace `Global.Q.Value`/`Global.PG.Q.Value` (library from the first MBR pass). |
| Points per peak `>= 6` across FWHM | community rule of thumb | Below this, quant precision and peak detection degrade; drives window/cycle design. |
| Fixed mass accuracy (MS2/MS1 ppm): timsTOF 15/15, Orbitrap Astral 10/4, TripleTOF 6600/ZenoTOF 20/20 | DIA-NN docs | Auto (0) is optimised on the first run and reused, so results depend on run order; fixed per-instrument values are preferred for production analyses. |
| Missed cleavages `1` | -- | Trypsin standard; raising it expands the search space and the multiple-testing burden. |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| `FileNotFoundError: report.tsv` or stale data | DIA-NN 1.9+ default is `report.parquet`, not `report.tsv` | Read `report.parquet` (`pd.read_parquet`); request legacy TSV explicitly only if needed |
| Matrix protein count < report count, looks like data loss | Matrices apply an extra 5% run-specific PG filter | Expected; do not compare the two counts as if equal |
| `-Inf` after log2 of the matrix | DIA-NN writes 0 for not-quantified | Convert `0 -> NaN` BEFORE log2/normalization |
| Cohort "hits" do not reproduce | Filtered per-run `Q.Value` only, not global | Filter `Global.PG.Q.Value <= 0.01` for the matrix |
| `easypqp convert --format diann` errors | No such option; `convert`/`library`/`insilico-library` are the real subcommands | `easypqp convert` per run, then `easypqp library ... *.psmpkl *.peakpkl`, or let FragPipe emit a DIA-NN-format library |
| `easypqp library`: `No PSMs files present` or `There is a psm.tsv but no peptide.tsv` | positional `*.psmpkl`/`*.peakpkl` from `easypqp convert` missing; `--psmtsv` given without `--peptidetsv` | run `convert` first; pass both tsvs and the pickles |
| `KeyError: 'report.pr.matrix.tsv'` | Matrix filename dotting varies by version (`pr_matrix` vs `pr.matrix`) | `ls` the output dir after a run and match the installed version's exact names |
| Noisy results on overlapping-window data in an engine without native support | Staggered data not demultiplexed | Demultiplex at conversion (peak picking first); DIA-NN supports overlapping windows natively |

## References

- Gillet LC, Navarro P, Tate S, et al. Targeted data extraction of the MS/MS spectra generated by data-independent acquisition: a new concept for consistent and accurate proteome analysis. *Mol Cell Proteomics* 2012;11(6):O111.016717.
- Rost HL, Rosenberger G, Navarro P, et al. OpenSWATH enables automated, targeted analysis of data-independent acquisition MS data. *Nat Biotechnol* 2014;32(3):219-223.
- Rosenberger G, Bludau I, Schmitt U, et al. Statistical control of peptide and protein error rates in large-scale targeted data-independent acquisition analyses. *Nat Methods* 2017;14(9):921-927.
- Searle BC, Pino LK, Egertson JD, et al. Chromatogram libraries improve peptide detection and quantification by data independent acquisition mass spectrometry. *Nat Commun* 2018;9:5128.
- Demichev V, Messner CB, Vernardis SI, Lilley KS, Ralser M. DIA-NN: neural networks and interference correction enable deep proteome coverage in high throughput. *Nat Methods* 2020;17(1):41-44.
- Meier F, Brunner AD, Frank M, et al. diaPASEF: parallel accumulation-serial fragmentation combined with data-independent acquisition. *Nat Methods* 2020;17(12):1229-1236.
- Amodei D, Egertson J, MacLean BX, et al. Improving precursor selectivity in data-independent acquisition using overlapping windows. *J Am Soc Mass Spectrom* 2019;30(4):669-684.
- Savitski MM, Wilhelm M, Hahne H, Kuster B, Bantscheff M. A scalable approach for protein false discovery rate estimation in large proteomic data sets. *Mol Cell Proteomics* 2015;14(9):2394-2404.
- Guzman UH, Martinez-Val A, Olsen JV, et al. Ultra-fast label-free quantification and comprehensive proteome coverage with narrow-window data-independent acquisition. *Nat Biotechnol* 2024;42:1855-1866.

## Related Skills

- spectral-libraries - Build experimental, chromatogram, or predicted libraries to search against
- quantification - Normalization, MaxLFQ roll-up, and matrix summarization after filtering
- differential-abundance - Moderated statistical testing of the protein matrix
- proteomics-qc - Per-run ID counts, missing-value rates, and acquisition QC
- data-import - Convert and load raw vendor/mzML data before the search
- peptide-identification - DDA spectrum-to-peptide matching (the non-DIA counterpart)
- workflows/proteomics-pipeline - End-to-end DIA-to-differential-abundance orchestration
