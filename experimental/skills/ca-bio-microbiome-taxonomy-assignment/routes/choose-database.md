# Which database and classifier

No script.

| Scenario | Use |
|----------|-----|
| 16S V4 (515F/806R), standard | SILVA region-matched (515-806) pre-trained naive Bayes, or DADA2 `assignTaxonomy` on SILVA |
| Primers known, no pre-trained artifact for them | `extract-reads` then `fit-classifier-naive-bayes` (RESCRIPt builds the reference) |
| Known habitat, best species accuracy | Weighted (clawback) naive-Bayes classifier |
| scikit-learn version error, or no retraining wanted | `classify-consensus-vsearch` (`routes/qiime2-vsearch.md`) |
| Conservative, novelty-aware, fewest over-calls | IDTAXA (`routes/idtaxa.md`) |
| One tree for 16S and shotgun | Greengenes2 (q2-greengenes2), closed-reference, GTDB-harmonized |
| Environmental or under-named bacteria | GTDB SSU reference via RESCRIPt (rank-normalized) |
| Fungal ITS | UNITE (species hypotheses) with naive Bayes or vsearch; never fixed-length-trim ITS |
| Protist or eukaryote 18S | PR2 or SILVA 18S |
| Reproducing an older study | RDP reference database |
| Species name from one 16S region | Do not; report genus, or `addSpecies` exact match only |
| Raw shotgun reads, not ASVs | `routes/other-skills.md` |

- Pick one database and release matched to the marker (16S: SILVA, GTDB or Greengenes2; ITS: UNITE; 18S: PR2 or SILVA). GTDB has no Eukarya, so it cannot classify 18S or ITS; a bacteria-only database gives empty or garbled ITS calls.
- Naive Bayes and vsearch-consensus perform comparably; confirm current best practice in the q2-feature-classifier documentation rather than hard-coding one.

Done when one database release and one classifier are named for the marker and region.
