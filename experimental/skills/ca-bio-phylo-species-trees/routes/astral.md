# Infer from existing gene trees
~~~bash
wastral -i gene_trees.nwk -o species_wastral.tre
~~~
Input is one Newick gene tree per line. Verify taxon identities and orthology; use an explicit species map for aliases or multiple individuals. Source recommends weighting or contracting weak standard-bootstrap branches (10% is a heuristic; UFBoot uses a different scale). For ASTER annotations:
~~~bash
astral -t 8 -u 2 -i gene_trees.nwk -o species_annot.tre
~~~
ASTER uses -t for threads and -u for annotation; classic Java ASTRAL uses different flags. Confirm executable version and reported species count. ASTER default branch units differ from wASTRAL; neither is time.
Done when mappings, support, units and input-tree uncertainty are reviewed.
