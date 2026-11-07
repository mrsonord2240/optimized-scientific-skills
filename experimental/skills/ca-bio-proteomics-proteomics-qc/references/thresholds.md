# Quantitative thresholds

Seed values for a per-instrument control chart, not standards-body limits. Status says how firm each is.

| Threshold | Status | Note |
|-----------|--------|------|
| Mass accuracy (internal/lock) median absolute error < 1-3 ppm, single mode centered on 0 | Convention | Width approaching the MS1 tolerance loses real IDs |
| iRT RT-fit R^2 > 0.99 (warn below) | Convention | Residual is LC noise on a stable gradient |
| FWHM alarm on a 20-30% rise over baseline; at least 8-10 points across the peak (floor about 5) | Published | Peak capacity tracks peptide IDs |
| % MS2 identified: < 20% bad, 20-35% ok, 35% and above great | PTXQC defaults | Generic marks, not a biological ceiling |
| Missed cleavages: 75-85% or more at 0; flag > 25-30% with 1 or more | Operational | Porcine trypsin is about 78% efficient even ideally |
| Replicate Pearson r (log2): technical > 0.98, biological 0.90-0.98, floor 0.8 | Convention | |
| Median CV (linear): technical < 10-20%, biological 20-40% | Convention | DIA below DDA |
| Completeness: valid in 50-70% or more of replicates in at least one condition | Convention | Filter before imputing |
| Perseus MNAR imputation: down-shift 1.8 SD, width 0.3 | Published | |
| Loading failure: total signal at or below 0.5x group median, or IDs below 0.8x | Script default | Judge on total signal and ID count |
| TMT channel: investigate beyond about 2x plex median, flag beyond 3-4x | Convention | |
| TMT isolation interference < 50% (< 30% stricter) | Convention | |
| DIA precursor and protein q both <= 0.01, global q, picked protein estimator | Standard | |
| Levey-Jennings: +/-2 SD warn, +/-3 SD action; QC every 4th-5th injection | Published | |
