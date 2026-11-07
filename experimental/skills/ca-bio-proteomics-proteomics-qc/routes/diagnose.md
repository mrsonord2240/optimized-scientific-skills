# Protein count dropped, cause unknown; queue drift

No script: this route localizes a fault from co-moving metrics. Read each metric with its co-readouts, never alone, and descend from the matrix to the instrument.

| Protein count drop with | Points to |
|-------------------------|-----------|
| Erratic TIC and injection time at its maximum | Spray or emitter |
| Retention time shifted, broad peaks (FWHM up), rising backpressure | Column |
| High contaminant fraction, normal spray and RT | Sample |
| High-charge (3+/4+) fraction up with missed cleavages up | Under-digestion |
| 1+ fraction up with short peptides | Ionization or spray, or contaminants |
| Bimodal mass-error histogram, normal ID rate | Monoisotopic mis-assignment or co-isolation: fix isotope-error tolerance and deisotoping, not calibration |

Level 1 evidence (instrument): median absolute mass error under 1-3 ppm, one peak centered on 0; iRT fit R^2 above 0.99; FWHM alarm on a 20-30% rise over baseline with at least 8-10 points across the peak. Tools: RawTools, RawBeans, rawrr, rawDiag, QuaMeter. DIA-NN reports go to `routes/diann-level1.md`.

- A fully tryptic peptide is mostly 2+; excess 3+/4+ comes from internal basic residues left by missed cleavages.
- Use interspersed QC injections every 4th-5th run and a Levey-Jennings chart: warn at +/-2 SD, act at +/-3 SD of the lab rolling baseline.
- Contaminants left in before normalizing shift the median; when their load differs by group the difference is normalized into real proteins. Strip them first (`routes/raw-signal.md`).

Done when the likely level (instrument, identification, sample) is named with the co-readouts that support it. Thresholds and their status: `references/thresholds.md`.
