## DECIPHER IDTAXA

**Goal:** Get a conservative, novelty-aware classification that refuses to descend into a clade the query likely does not belong to.

**Approach:** Convert ASV sequences to a DNAStringSet, classify with a pre-trained DECIPHER trainingSet, then flatten the per-rank output to a matrix, mapping IDTAXA's "unclassified_" placeholders to NA.

```bash
# Pre-trained trainingSet (.RData providing `trainingSet`, or .rds):
Rscript scripts/idtaxa_classify.R seqtab_nochim.rds SILVA_SSU_r138_2019.RData taxa_idtaxa.rds
# No pre-trained set for the marker/region: train with LearnTaxa() from a reference FASTA + matching
# "Root;domain;phylum;...;genus;" taxonomy strings (one per sequence, same order), then classify:
Rscript scripts/idtaxa_classify.R seqtab_nochim.rds trainingset.rds taxa_idtaxa.rds region-matched-ref.fasta region-matched-ref-taxonomy.txt
# Append `rerun` to repeat the seeded IdTaxa() call once and require identical() output.
```

`scripts/idtaxa_classify.R` (checked on DECIPHER 3.2.0) writes an ASV x rank character matrix (domain..species, NA where unclassified). What it does, and why:

- **Seed before EVERY `LearnTaxa()` and `IdTaxa()` call** (`set.seed(100)`; any fixed integer works). Both are internally stochastic: unseeded, two `LearnTaxa()` calls on identical input give non-identical trainingSets, and two `IdTaxa()` calls on the identical trainingSet and query set differ at ~3-4% of genus calls (a larger margin than `assignTaxonomy()`). Seeded, `IdTaxa()` is bit-identical at every rank including genus, with `processors=NULL` (multithreaded) and `processors=1`, and `LearnTaxa()` gives an `identical()` trainingSet.
- `threshold = 60` is the DECIPHER default confidence cutoff; raise it for stricter calls. IDTAXA stops descending (leaves the rank unclassified) when the query likely belongs to a taxon absent from the reference - the intended anti-over-classification behaviour.
- **Flatten positionally, never by name.** `x$rank` is populated only when the trainingSet was built with `LearnTaxa()`'s optional `rank=` data.frame (a 5-column Index/Name/Parent/Level/Rank table, rarely available outside DECIPHER's pre-built sets). Without it `x$rank` is NULL and `match(ranks, x$rank)` silently returns all-NA at every rank with no error. `x$taxon[1]` is always "Root" and the rest are domain..genus/species in order, so the script drops "Root", maps `unclassified_*` placeholders to NA and pads to the fixed rank depth.
- **Sanity check:** the script stops if every genus call is NA (a wrong-region/wrong-marker trainingSet or a flattening bug gives all-NA with no error).
- **MEMORY:** `LearnTaxa()` against a full, un-subsampled reference (400K+ sequences) needs tens of GB of RAM and can crash on constrained hardware; subsample the reference (e.g. ~60,000 sequences) if it does.
