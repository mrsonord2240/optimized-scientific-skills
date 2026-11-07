# Draw a simple tree
~~~bash
python scripts/draw.py tree.nwk --ascii
python scripts/draw.py tree.nwk --output tree.pdf --title "Phylogram" --support-label "SH-aLRT / UFBoot" --ladderize
~~~
Done when the ASCII view or new vector figure is inspected, branch-length units and support measure are named, and layout/root assumptions are documented. Supply the actual support measure rather than copying the example blindly. Ladderization changes display order, not biology. The script refuses over 200 tips; Bio.Phylo supports rectangular layouts only. Choose SVG/PDF/EPS for vector output; if raster is required use at least 600 dpi through an adapted source recipe. See [guidance](../references/guidance.md).
