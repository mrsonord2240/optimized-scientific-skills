# Filter host organelle and off-target features

```bash
bash scripts/filter_organelle.sh table.qza taxonomy.qza rep-seqs.qza filtered_out
```

Run after assignment and before diversity or differential abundance. Universal 16S primers amplify host mitochondrial and plastid 16S; left in, they inflate the table and deflate every real taxon by compositional closure. The script runs `qiime taxa filter-table` and `filter-seqs` with `--p-exclude mitochondria,chloroplast`, a case-insensitive substring match on the lineage that works for any database. Heaviest in plant, rhizosphere and host-tissue or biopsy samples (often most reads): inspect per-sample read retention afterwards.

phyloseq form, SILVA 138 only (Chloroplast at Order, Mitochondria at Family); for GTDB, Greengenes2 or RDP check the rank first with `get_taxa_unique(ps, 'Order')`:

```r
ps <- subset_taxa(ps, is.na(Order)  | Order  != 'Chloroplast')
ps <- subset_taxa(ps, is.na(Family) | Family != 'Mitochondria')
```

Exception: in phototroph, aquatic or microbial-mat communities the Chloroplast bin can be the signal (cyanobacterial and algal-plastid 16S are hard to separate). Inspect what falls in it before excluding.

Done when the filtered table and sequences exist and per-sample retention was checked.
