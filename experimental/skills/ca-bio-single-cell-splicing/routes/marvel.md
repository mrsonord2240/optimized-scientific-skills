# MARVEL plate data

```bash
Rscript scripts/marvel_psi.R annotation.gtf rmats_out star_pass2 cells.rds intron_counts.tsv gene_features.tsv tpm.tsv marvel.rds SE
```

The inputs and filenames follow `references/marvel.md`; all five rMATS event tables must exist. The wrapper saves an object with PSI for the requested event class using CoverageThreshold 10. For RI append the actual read length after RI; it uses two threads. Match contig names and check coverage and alignment. Modality and group-comparison code remains input-dependent in the reference; adapt it to real cell labels and interpret comparisons with donor replication.

Done when the requested PSI object and alignment checks are available, and any requested comparisons are reviewed.
