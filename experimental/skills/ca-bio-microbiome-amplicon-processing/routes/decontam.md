# Remove reagent contaminants (decontam)

```bash
Rscript scripts/decontam_filter.R dada2_out/seqtab_nochim.rds metadata.csv decontam_out [threshold]
```

`metadata.csv`: first column is the sample name; `is_control` (TRUE for extraction blanks and no-template PCR) and/or `dna_conc` (qPCR or Qubit). Method follows the columns: both gives combined, only controls gives prevalence, only concentration gives frequency. Threshold 0.1 by default; 0.5 is aggressive. Writes `contaminants.csv` and `seqtab_clean.rds` and prints the ASVs and reads removed; report that.

- Low-biomass samples (skin, biopsy, BAL, sterile-site swabs) can be dominated by the kit contaminants. Never interpret one without negative controls and, ideally, a positive mock community. A plausible "community" in a near-sterile sample, reagent-associated genera, or results that track DNA yield are the signature.
- Do not substitute a flat relative-abundance cutoff (drop ASVs below 1% or 5%). Contaminant status comes from controls and DNA concentration, not abundance: in a synthetic test the kit contaminant (*Ralstonia*, 5.1% of reads) outnumbered a real member (*Escherichia-Shigella*, 4.3%), so any cutoff removing one removes the other or keeps both.
- The table must be samples in rows, ASVs in columns, as `dada2_asv.R` writes it.

Done when `seqtab_clean.rds` exists and the removed ASV and read counts are stated.
