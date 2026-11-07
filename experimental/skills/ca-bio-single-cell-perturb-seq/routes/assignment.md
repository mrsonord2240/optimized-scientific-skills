# Guide-count mixture assignment
~~~sh
python scripts/assign_guides.py guide_counts.h5ad assigned.h5ad
~~~
Input is cells by guides with raw counts in X. The source mixture method requires its optional JAX backend and writes assigned_guide metadata; check package version and actual outputs. Inspect ambient contamination, controls, doublets and multi-guide calls. A threshold fallback is a distinct method and should not be silently substituted.
**Done when:** assignment, controls, uncertainty and per-guide distributions are reviewed and saved.
