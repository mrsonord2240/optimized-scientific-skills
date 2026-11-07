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
| Discovery cohort, no wet-lab library | Generate a predicted library from FASTA **without raw files**, then search against it | Bounded, better-calibrated search vs raw directDIA; the modern default. Required by DIA-NN 2.x. |
| Quick single-run discovery, Astral 2-Th data | DIA-NN library-free (directDIA) | Near-non-chimeric MS2 makes directDIA trustworthy |
| Have a deep experimental/chromatogram library | `diann --lib library.parquet --fasta db.fasta` (no `--fasta-search`) | Targeted extraction is most sensitive when the library matches |
| Need auditable run/experiment/global FDR for a regulated submission | OpenSWATH + PyProphet | Explicit q-value contexts per Rosenberger 2017 |
| Building chromatogram-library depth for one project on Orbitrap | EncyclopeDIA (GPF) -> .elib -> DIA-NN | Empirical RT and real fragmentation in the project's own LC |
| Large cohort (hundreds-thousands of runs) | Empirical library from 20-100 high-quality runs (library-free search + `--gen-spec-lib`), then search ALL runs with that library, MBR off, fixed mass accuracies; filter on `Global.Q.Value`/`Global.PG.Q.Value` | DIA-NN docs: MBR is a convenience feature (second pass held in memory); one fixed library and fixed tolerances keep batches comparable |
| Staggered/overlapping acquisition | DIA-NN directly (native support); for engines without it, demultiplex at conversion with peak picking first | Demultiplexing is an optional few-% gain in DIA-NN; required where the engine assumes non-overlapping windows |
| PTM / peptidoform-resolved work | DIA-NN `--peptidoforms` + matched variable mods | Peptidoform-resolved target-decoy scoring |

Default when uncertain: use the two-stage DIA-NN predicted-library route in [`examples/diann_analysis.sh`](examples/diann_analysis.sh): generate the predicted library from FASTA with **no raw files**, then search raw files against that library with FASTA digest and prediction off. Fix mass accuracies for the instrument, then filter `Q.Value`, `PG.Q.Value`, `Global.Q.Value` and `Global.PG.Q.Value` at 0.01 for cross-run matrices (DIA-NN 1.9.x with MBR: `Lib.Q.Value`/`Lib.PG.Q.Value` instead of the `Global.*` pair).


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

An EasyPQP `library.tsv` is not automatically a DIA-NN library. With DIA-NN 2.6.1, a TSV must include fragment-annotation fields `FragmentCharge`, `FragmentType`, and `FragmentSeriesNumber`; EasyPQP's generic TSV did not provide them in the audited route. Prefer FragPipe's DIA-NN-format export. If converting another TSV, validate those fields and a small DIA-NN load before a cohort search; do not add them by guessing fragment annotations.

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
| `WARNING: incorrect settings, the in silico-predicted library must be generated in a separate pipeline step` | DIA-NN 2.x was given raw files together with FASTA digest / prediction | Generate the `.predicted.speclib` from FASTA with no raw files, then run a separate `--lib` search with `--fasta-search` and `--predictor` off |
| DIA-NN cannot load an EasyPQP `library.tsv`, or reports missing fragment fields | The TSV lacks DIA-NN fragment annotations (`FragmentCharge`, `FragmentType`, `FragmentSeriesNumber`) | Prefer FragPipe's DIA-NN-format export; otherwise use a validated converter and test-load the library before searching the cohort |
| Noisy results on overlapping-window data in an engine without native support | Staggered data not demultiplexed | Demultiplex at conversion (peak picking first); DIA-NN supports overlapping windows natively |
