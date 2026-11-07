# Errors and wrong-looking results

| Symptom | Cause | Fix |
|---------|-------|-----|
| "The scikit-learn version (X) used to generate this artifact does not match ... (Y)" | Pre-trained `.qza` built under another QIIME2 release; the classifier is a pickled scikit-learn model | Download the classifier built for the installed release, retrain with `fit-classifier-naive-bayes`, or use `routes/qiime2-vsearch.md`. Never pin a hard-coded classifier URL |
| Systematically shallow or wrong calls versus a region-matched run | Full-length-trained classifier on a sub-region such as V4 | Use a primer-matched classifier or `extract-reads` then train |
| Confident species labels that `addSpecies` would leave NA | Over-reading species from one 16S region | Report genus; species via exact match or IDTAXA only |
| Spurious genus mismatches across cohorts; empty ITS calls | Different databases (SILVA, GTDB), or a bacteria-only database for ITS or 18S | One database and release matched to the marker (`routes/choose-database.md`) |
| Over-deep labels, or silently dropped or force-filled Unassigned | Confidence default left unexamined | State the threshold; keep the truncated output (`routes/rank-and-confidence.md`) |
| Large read fraction labelled Mitochondria or Chloroplast | Host organelle 16S amplified by universal primers | `routes/organelle.md` |
| OOM kill during `classify-sklearn` | Large classifier times `--p-n-jobs` | Lower `--p-n-jobs` or `--p-reads-per-batch`; use a region-extracted, smaller classifier |
| All Unassigned at domain level | Off-target ASVs (host, chimera, primer artifact) or wrong-orientation reads | Filter off-target; leave read orientation on `auto`; document, do not force-fill |
| Empty or garbage ITS taxonomy | GTDB, or position-trimmed ITS | UNITE; remove ITS primers but never fixed-length-trim ITS |
| IDTAXA flattening returns all-NA at every rank, no error | `x$rank` is NULL; flattened by name | Flatten positionally (the script does) |
| Two runs of `assignTaxonomy()`, `IdTaxa()` or `LearnTaxa()` differ on identical input | Stochastic, no seed | `set.seed()` before every call |

Done when the failing command has run to completion and the cause is addressed.
