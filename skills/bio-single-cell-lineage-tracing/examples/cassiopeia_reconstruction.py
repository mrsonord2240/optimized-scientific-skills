'''Lineage tree reconstruction with Cassiopeia'''
# Reference: cassiopeia 2.0.0, numpy 1.26.4 | Verify API if version differs
# Run in a dedicated Cassiopeia env (numpy<2) -- see SKILL.md "Installation and
# Version Compatibility" for why Cassiopeia and CoSpar/scanpy cannot share one env.
import cassiopeia as cas
import pandas as pd
import numpy as np

# Option 1: Load pre-built character matrix
# Rows = cells, Columns = barcode sites
# Values: 0 = unedited, 1-N = different mutations, -1 = missing
char_matrix = pd.read_csv('character_matrix.csv', index_col=0)
cell_meta = pd.read_csv('cell_metadata.csv', index_col=0)

# Create CassiopeiaTree object
tree = cas.data.CassiopeiaTree(
    character_matrix=char_matrix,
    cell_meta=cell_meta
)

print(f'cells {tree.n_cell}  characters {tree.n_character}  missing {(char_matrix.values == -1).mean():.2%}')

# HybridSolver is the scalable default: greedy top split, exact solve on small subclades.
# bottom_solver=ILPSolver() needs a separately licensed Gurobi install (gurobipy); this
# uses VanillaGreedySolver() as a license-free fallback -- swap in ILPSolver() if you
# have a Gurobi license, for a near-optimal solve on the small subclades instead.
# collapse_mutationless_edges removes internal edges with no supporting mutation
solver = cas.solver.HybridSolver(
    top_solver=cas.solver.VanillaGreedySolver(),
    bottom_solver=cas.solver.VanillaGreedySolver(),
    cell_cutoff=200
)
solver.solve(tree, collapse_mutationless_edges=True)

newick = tree.get_newick()

# Infer ancestral states on internal nodes
tree.reconstruct_ancestral_characters()

# Compare against a neighbor-joining tree to gauge topology robustness
# add_root=True is required by the installed API -- omitting it raises DistanceSolverError
nj_tree = cas.data.CassiopeiaTree(character_matrix=char_matrix, cell_meta=cell_meta)
cas.solver.NeighborJoiningSolver(
    dissimilarity_function=cas.solver.dissimilarity_functions.weighted_hamming_distance,
    add_root=True
).solve(nj_tree)
rf, rf_max = cas.critique.robinson_foulds(tree, nj_tree)
print(f'Robinson-Foulds {rf}/{rf_max}  triplets-correct {cas.critique.triplets_correct(tree, nj_tree)}')

# Installed Cassiopeia 2.0.0's cas.pl only exposes iTOL export (upload_and_export_itol,
# needs an iTOL API key), not a local matplotlib plot function -- plot `newick` with an
# external tool (e.g. ete3, already a Cassiopeia dependency, or iTOL) if you need a figure.
