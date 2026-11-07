# DECIPHER IDTAXA

```bash
# pre-trained trainingSet (.RData providing `trainingSet`, or .rds)
Rscript scripts/idtaxa_classify.R seqtab_nochim.rds SILVA_SSU_r138_2019.RData taxa_idtaxa.rds
# no pre-trained set for the marker or region: train with LearnTaxa(), then classify
Rscript scripts/idtaxa_classify.R seqtab_nochim.rds trainingset.rds taxa_idtaxa.rds region-matched-ref.fasta region-matched-ref-taxonomy.txt
```

Training taxonomy: one `Root;domain;phylum;...;genus;` string per reference sequence, same order. Append `rerun` to repeat the seeded `IdTaxa()` once and require `identical()` output. The script (checked on DECIPHER 3.2.0) writes an ASV x rank character matrix (domain to species, NA where unclassified).

- Seed before every `LearnTaxa()` and `IdTaxa()` call (`set.seed(100)`; any fixed integer). Both are stochastic: unseeded, two `IdTaxa()` calls differ at about 3-4% of genus calls and two `LearnTaxa()` calls give different trainingSets. Seeded, results are identical at every rank, with `processors=NULL` or `1`.
- `threshold = 60` is the default confidence cutoff; raise it for stricter calls. IDTAXA stops descending, leaving the rank unclassified, when the query probably belongs to a taxon absent from the reference.
- Flatten positionally, never by name. `x$rank` is NULL unless the trainingSet was built with `LearnTaxa(rank=)`, and `match(ranks, x$rank)` then returns all-NA with no error. The script drops `Root`, maps `unclassified_*` to NA and pads to seven ranks.
- The script stops if every genus call is NA (wrong region or marker, or a flattening bug).
- `LearnTaxa()` on a full reference (400K+ sequences) needs tens of GB of RAM; subsample (about 60,000 sequences) if it crashes.

Done when the matrix is saved, genus calls exist, and the trainingSet's database release and region are stated.
