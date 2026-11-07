# How deep a rank, and which threshold

No script.

- A single hypervariable region (V4, about 250 bp) is identical across many species and often genera; the discriminating substitutions lie elsewhere in the gene. Classifiers still emit a species name. Report the rank the data supports: genus at best for many bacteria, family for poorly resolved clades. Mostly genus or family with few species is correct behaviour; do not lower confidence to force species.
- ITS resolves to species far more reliably than 16S. 18S resolves coarsely, like 16S.
- Confidence is certainty given that the taxon is in the reference. If the true organism is absent, the classifier confidently returns the nearest wrong relative.
- A threshold trades sensitivity for specificity: lowering over-classifies (deeper but wronger); raising truncates to shallower, more reliable ranks. State it, tune it to region and database, and keep the unassigned output as is.

| Setting | Default | Meaning |
|---------|---------|---------|
| `classify-sklearn --p-confidence` | 0.7 | Below it the lineage is truncated to a more confident rank |
| `assignTaxonomy` `minBoot` | 50 (80 common) | Ranks below the floor are NA |
| IDTAXA `threshold` | 60 | Per-node confidence cutoff; raise for stricter calls |
| `classify-consensus-vsearch` `--p-min-consensus` | 0.51 | Majority of accepted hits must agree on a rank |
| `classify-consensus-vsearch` `--p-perc-identity` | 0.8 | Minimum query-reference identity for a hit |
| `classify-consensus-vsearch` `--p-maxaccepts` | 10 | Top hits kept per query |
| `addSpecies` | exact (100%) | Species only from an exact amplicon match; high precision, low recall |
| `extract-reads --p-min-length` | 50 | Drops too-short in-silico amplicons that would mistrain the classifier |

Done when each reported label's rank is justified by the region and the threshold is stated.
