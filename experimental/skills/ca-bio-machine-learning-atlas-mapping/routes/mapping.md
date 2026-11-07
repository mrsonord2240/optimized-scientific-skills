# scANVI reference and query mapping
~~~sh
python scripts/scarches_annotation.py reference.h5ad query.h5ad annotated.h5ad
~~~
Both objects need raw counts in layers/counts and batch metadata; reference also needs cell_type labels and aligned gene IDs. This source template trains a reference scVI/scANVI model and fine-tunes a query, using seed 0, 100/20/100 epochs and weight_decay 0.0. It is not a saved-model-only projection. Review reference coverage, QC, feature overlap, training behavior and package APIs before use. The 0.2 weighted-neighbor disagreement gate is a source example, not a universal out-of-distribution detector.
**Done when:** output, reference provenance, settings, uncertainty and independent label checks are reviewed.
