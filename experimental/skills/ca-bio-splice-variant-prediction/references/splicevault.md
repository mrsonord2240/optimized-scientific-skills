## SpliceVault for Empirical Mis-Splicing Outcomes

**Goal:** Predict the *type* of mis-splicing (exon skipping vs cryptic site activation) given a canonical-disrupting variant.

**Approach:** Query the public SpliceVault table (GRCh38 SNVs at annotated splice sites, one row per transcript) by remote tabix; `scripts/splicevault_lookup.py` returns the parsed Top-4 events. Needs `pysam` and `pandas`; downloading the 884 MB `SpliceVault_data_GRCh38.tsv.gz` is optional (`--tsv local.gz`). Bulk SQL dumps at `storage.googleapis.com/misspl-db-data` are requester-pays (billing project required); the web portal is `kidsneuro.shinyapps.io/splicevault`; the same table backs the SpliceVault VEP plugin (Ensembl 111+).

```bash
curl -O https://ftp.ensembl.org/pub/current_variation/SpliceVault/SpliceVault_data_GRCh38.tsv.gz.tbi   # index only; rows stream over HTTP
python splicevault_lookup.py chr17 7674292 T C --transcript ENST00000269305 --tbi SpliceVault_data_GRCh38.tsv.gz.tbi
# ENST00000269305 Acceptor_loss at chr17:7674291; SpliceAI delta 1; out-of-frame Frameshift:3/4; 199336 samples
# Top1 CA +47 0.4% Frameshift | Top2 CA -50 0.08% Frameshift | Top3 ES 7 0.03% Frameshift | Top4 CA -70 0.03% inFrame
```

Event types: ES = exon skipping (impact = skipped exon number), CA/CD = cryptic acceptor/donor (impact = nt from the annotated site, transcript sense). Checked against SpliceAI on the same variant: Top1 `CA +47` is SpliceAI's acceptor gain (DS_AG 0.85, DP_AG +47). `out_of_frame` is the fraction of the Top-N events that shift the frame; >= 3/4 in-frame suggests a variant that may be splice-rescued rather than LoF. An empty result means the variant is not a catalogued splice-site SNV (or REF/ALT/build is wrong).

SpliceVault (Dawes 2023 *Nat Genet*; 335,663 RNA-seq samples) reports that the **300K-RNA Top-4 events** predict variant-associated mis-splicing with 92% sensitivity, identifying 96% of exon-skipping events and 86% of cryptic-site events in 140 clinical cases with RNA testing. Use it when the question is not "will splicing change?" but "what specific aberrant splicing will occur?".
