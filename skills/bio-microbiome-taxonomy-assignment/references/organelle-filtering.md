## Filtering Host Organelle and Off-Target Features

**Goal:** Remove host mitochondrial 16S, chloroplast/plastid 16S, and domain-unassigned features BEFORE diversity and differential abundance - universal 16S primers amplify host organelle rRNA, and leaving it in inflates the feature table and deflates every real taxon by compositional closure.

**Approach:** Use the taxonomy just assigned to exclude the Mitochondria and Chloroplast lineages (the labels enable the filter), then carry the filtered table and sequences forward.

```bash
# QIIME2: exclude mitochondria + chloroplast (case-insensitive substring match on the lineage)
qiime taxa filter-table --i-table table.qza --i-taxonomy taxonomy.qza \
    --p-exclude mitochondria,chloroplast \
    --o-filtered-table table-no-organelle.qza
qiime taxa filter-seqs --i-data rep-seqs.qza --i-taxonomy taxonomy.qza \
    --p-exclude mitochondria,chloroplast \
    --o-filtered-data rep-seqs-no-organelle.qza
```

```r
# phyloseq (SILVA ranks: Order 'Chloroplast', Family 'Mitochondria'); the is.na guard keeps unranked taxa
ps <- subset_taxa(ps, is.na(Order)  | Order  != 'Chloroplast')
ps <- subset_taxa(ps, is.na(Family) | Family != 'Mitochondria')
```

Organelle contamination is heaviest in plant, rhizosphere, and host-tissue/biopsy samples (often the majority of reads); inspect per-sample read retention after filtering. Exception: in phototroph-focused, aquatic, or microbial-mat communities, Chloroplast-binned 16S can be the signal of interest (cyanobacterial vs algal-plastid 16S are hard to separate) - inspect what falls in the Chloroplast bin before excluding it. The phyloseq rank-equality form is SILVA-138-specific (Chloroplast at Order, Mitochondria at Family); for GTDB/Greengenes2/RDP verify the rank (`get_taxa_unique(ps, 'Order')`) or use the QIIME2 substring exclude, which is reference-robust.
