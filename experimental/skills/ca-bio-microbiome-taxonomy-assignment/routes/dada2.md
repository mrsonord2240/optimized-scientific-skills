# DADA2 assignTaxonomy + addSpecies

```bash
Rscript scripts/dada2_assign.R seqtab_nochim.rds silva_nr99_v138.1_train_set.fa.gz taxa.rds silva_species_assignment_v138.1.fa.gz
```

Arguments: ASV table (ASV sequences as columns), DADA2-formatted training FASTA, output `.rds`, optional species FASTA, optional trailing `rerun`. `MIN_BOOT` env sets `minBoot` (default 50). The script seeds, runs `assignTaxonomy(minBoot, tryRC=TRUE)` to genus, stops if no ASV gets a genus, runs `addSpecies` when a species FASTA is given, prints the assigned fraction per rank, and records the training file and `minBoot` as attributes of the saved matrix. ITS: use a UNITE DADA2 reference and never position-trim ITS.

- `minBoot` 50 is the DADA2 default and the RDP recommendation for reads up to 250 nt; 80 is a stricter choice. Ranks below it come back NA.
- `addSpecies` assigns by exact (100%) match only. A species is reported when the ASV is identical to a reference over the amplicon, else NA.
- All-NA output with no error means a wrong-region or wrong-marker reference, bad orientation, or mis-formatted FASTA; the script stops on zero genus calls.
- Seeded runs agree at genus; about 2 of 770 cells at Kingdom or Order can still differ, with `multithread=FALSE` too. Add `rerun` to repeat once and require identical genus calls.
- The reference should be trimmed to the amplicon region; a full-length SILVA training set on V4 reads degrades calls (`routes/errors.md`).

Done when `taxa.rds` is written, genus is assigned for a plausible fraction of ASVs, and the report names the database release and `minBoot`.
