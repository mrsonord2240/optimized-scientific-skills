# Convert plain tree files
~~~bash
python scripts/convert_plain.py tree.nwk tree.xml --from-format newick --to-format phyloxml
~~~
Use only when Bio.Phylo's format behavior is suitable. Extract meaningful BEAST/MrBayes/NHX annotations first with [annotations](annotations.md); conversion can corrupt or drop metadata. Use distinct input/output paths and a fresh destination. Check all tip labels, tree counts, support fields, rootedness and required annotations after round-trip.
Done when the destination and its retained or lost fields are documented.
