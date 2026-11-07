# Pathogenicity predictors for PP3/BP4

No script. Use exactly one calibrated predictor per evidence type (missense; splicing) at its calibrated strength. Stacked predictors are not independent: REVEL is a random forest over 13 component scores including SIFT and PolyPhen, and most predictors train on ClinVar or HGMD labels, so benchmarks on those databases are partly circular.

| Predictor | Scope | Use |
|-----------|-------|-----|
| REVEL | rare missense | Best-calibrated single missense tool. PP3_Supporting at 0.644 or above; BP4_Supporting at 0.290 or below |
| AlphaMissense | missense, proteome-wide | Not trained on ClinVar labels. Use the current ClinGen SVI calibrated thresholds, not the developer class cutoffs |
| CADD | all variant types | Genome-wide and non-coding ranking, not missense PP3 |
| SIFT, PolyPhen-2 | missense | Not standalone evidence |
| SpliceAI | splicing | The splicing predictor |

- The REVEL Moderate and Strong cutoffs commonly quoted (PP3_Moderate 0.773 or above, PP3_Strong 0.932 or above; BP4_Moderate 0.183 or below, BP4_Strong 0.016 or below) are not uniformly reproduced. Verify against the Pejaver 2022 supplement or the current ClinGen SVI table before hard-coding them. PP3 and BP4 are mutually exclusive.
- At developer thresholds SIFT "deleterious" and PolyPhen "probably damaging" did not reach Supporting for PP3; calibrated PolyPhen-2 intervals do (Supporting [0.978, 0.999), Moderate 0.999 or above). Both call many rare missense damaging and are inside REVEL, so quoting them with it double-counts.
- CADD 20 or above is not "pathogenic": developer CADD 20 mapped to Moderate evidence for benign. Calibrated CADD reaches PP3_Supporting at [25.3, 28.1) and Moderate at 28.1 or above.
- SpliceAI gives four deltas (acceptor gain/loss, donor gain/loss), each 0-1; the maximum is the headline. Developer cutoffs: 0.2 high recall, 0.5 recommended, 0.8 high precision. Know whether the run is masked or raw; the default window is +/-50 bp, so deep-intronic effects need a wider window (up to +/-10 kb); a delta needs the ClinGen splicing calibration to become PS3 or PP3 strength. Pangolin is a tissue-aware alternative; check current ClinGen splicing guidance.
- The AlphaMissense and splicing calibrations are still evolving: verify the current ClinGen SVI approved list before standardizing.

Done when one predictor per evidence type is named with its version and the calibrated strength applied.
