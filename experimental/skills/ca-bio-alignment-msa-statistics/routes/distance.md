# Corrected distances and uncertainty

```bash
python scripts/statistics.py identity alignment.fasta identity.json --method pid4
```

For model-corrected distances read `references/distance-correction.md`; exploratory identity is not a corrected distance. For important phylogenetic or selection work, assess column confidence with the multiple-alignment skill before inference. Protein identity below 25%, many gappy columns, or regional identity shifts call for checking homology and alignment quality.

For the source exploratory Kimura protein correction run `python examples/kimura_protein_distance.py alignment.fasta`; values at or above 85% observed divergence are reported as saturated. For confidence and model selection use the canonical skills `bio-alignment-multiple` and `bio-phylo-modern-tree-inference`.

Done when the requested distance or uncertainty assessment is reported with its model and limitations.
