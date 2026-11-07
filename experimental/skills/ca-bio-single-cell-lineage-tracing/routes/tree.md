# Reconstruct a scar-character tree
~~~sh
python scripts/reconstruct_tree.py character_matrix.csv tree.newick
~~~
Rows are cells and columns are character sites: 0 is unedited, positive integers are states, and -1 is missing. Preserve that distinction and assess missingness, doublets, saturation, collisions and homoplasy first. The source uses VanillaGreedySolver and collapses mutationless edges; this does not weight indel priors. Compare independent solvers and quantify topology uncertainty before interpretation. Check installed Cassiopeia API.
**Done when:** Newick, input provenance, solver, missingness and robustness limitations are documented.
