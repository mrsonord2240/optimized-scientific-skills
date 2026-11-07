# Per-label marker agreement
~~~sh
python scripts/label_marker_check.py --reference reference.h5ad --query annotated.h5ad --markers curated_markers.json --out marker_check.tsv
~~~
Both objects need counts layers; reference has cell_type and query has predicted_label. Curated JSON maps labels to marker lists. The source defaults require 20 cells and flag retained detection below 0.6; calibrate for the tissue and gene set. A passed group-level check is not proof of identity, and unchecked labels are unverified. Read [failure modes](../references/failure-modes.md).
**Done when:** flags, unchecked labels and their counts are recorded with independent evidence.
