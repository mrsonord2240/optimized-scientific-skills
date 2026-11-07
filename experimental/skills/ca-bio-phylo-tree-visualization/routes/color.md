# Color a rooted clade
~~~bash
python scripts/draw.py rooted.nwk --output clade.pdf --clade Human Chimp Gorilla --title "Rooted phylogram"
~~~
Done when the intended tip set exactly matches the MRCA descendants and the output is visually checked. Root first using bio-phylo-tree-manipulation; the input's rooting is a scientific assumption, not verified by this script. The clade is red and inherited by descendants. State the branch-length units and verify support assignments after rerooting.
