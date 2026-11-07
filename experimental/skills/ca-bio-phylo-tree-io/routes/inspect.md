# Inspect tree files
~~~bash
python scripts/inspect_trees.py trees.nwk --format newick
~~~
Uses Bio.Phylo parse for one or more trees, prints tips, rooted flag and internal name/confidence/comment fields through UTF-8. Numeric or dual-support labels require source-program interpretation; comments are opaque, not structured BEAST metadata. Nexus translation and rootedness need verification.
Done when tree count, labels, support representation and rootedness are reviewed.
