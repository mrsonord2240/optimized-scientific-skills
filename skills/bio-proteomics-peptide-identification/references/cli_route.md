# Command-line DDA route and rescoring

Read when running Sage, Comet or MS-GF+ from the command line, rescoring the `.pin` with Percolator, or pooling runs. Decoy tags and the database build stay in `SKILL.md`; `examples/dda_search.sh` runs this route end to end.

### Run a DDA Search from the Command Line

**Goal:** Search a centroided mzML against the database and emit PSMs plus Percolator features.

**Approach:** All three engines below run one concatenated search and write a Percolator `.pin`, so the rescoring step is the same for each. Convert vendor raw first (`msconvert --mzML --zlib --filter "peakPicking vendor msLevel=1-"` -> data-import). Set tolerances to the instrument, not to the default: 10 ppm precursor and high-res fragment settings for an Orbitrap, 0.6 Da / ion-trap binning for CID. `examples/dda_search.sh` runs the whole route (database -> search -> Percolator -> 1% list) for `ENGINE=sage` or `ENGINE=comet`; `MZML` may list several runs, which are searched and rescored together, and it stops before Percolator if the pin's `Label` column holds no decoys (the decoy-tag trap in `SKILL.md`).

```bash
# Sage 0.14.6 -- generates its own 'rev_' decoys, so it takes the TARGET-ONLY FASTA.
# search.json: {"database": {"enzyme": {"missed_cleavages": 2, "min_len": 7, "max_len": 30,
#   "cleave_at": "KR", "restrict": "P"}, "static_mods": {"C": 57.0215},
#   "variable_mods": {"M": [15.9949]}, "max_variable_mods": 2,
#   "decoy_tag": "rev_", "generate_decoys": true, "fasta": "human.fasta"},
#  "precursor_tol": {"ppm": [-10, 10]}, "fragment_tol": {"ppm": [-20, 20]},
#  "isotope_errors": [0, 1], "deisotope": true, "predict_rt": true, "report_psms": 1}
sage search.json sample.mzML --write-pin     # -> results.sage.tsv + results.sage.pin

# Comet 2026.02 -- searches the concatenated DB built in SKILL.md.
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

Rescoring pays where the engine's own score is weakest. On one Orbitrap Astral 5-min DDA run (250 pg HYE load, 6,135 MS2 spectra, 31,437-protein database, PXD070049), Comet's raw `-log10(e-value)` gave **916** PSMs at 1% FDR and Percolator lifted the same search to **1,144** (+25%); Sage's `sage_discriminant_score` is already a learned score, so Percolator moved it from **1,406** to **1,398** (-0.6%) -- no gain to be had. Sage 0.14.6 wins this input outright; Comet 2026.02 with Percolator lands 19% behind it, and MS-GF+ v2024.03.26 read straight off `-log10(SpecEValue)` with no rescoring gives **656**, because a calibrated E-value buys cross-instrument comparability, not raw yield. **Do not generalise these counts**: one run, one low-load short-gradient method, each engine at its own idiomatic high-res settings. Rescoring also needs training data -- pooling the three DDA runs (6,864 PSMs) gave Percolator 5,006 PSMs against Sage's own 4,961, while on a single run of 1,939 PSMs it had too few positives to improve anything. To pool, give the engine every run (Sage takes several mzML in one call; for Comet concatenate the runs' pins with the header kept once) and rescore once.
