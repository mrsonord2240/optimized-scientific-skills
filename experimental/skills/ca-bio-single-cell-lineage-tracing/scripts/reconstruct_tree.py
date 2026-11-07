"""Reconstruct an illustrative Cassiopeia tree from an existing character matrix."""
import argparse
from pathlib import Path
import pandas as pd
import cassiopeia as cas

p = argparse.ArgumentParser(description=__doc__)
p.add_argument("matrix_csv")
p.add_argument("output_newick")
args = p.parse_args()
matrix = pd.read_csv(args.matrix_csv, index_col=0)
if not matrix.index.is_unique or matrix.empty or matrix.isna().any().any():
    p.error("matrix needs unique cell IDs and explicit -1 missing states")
if ((matrix < -1) | (matrix % 1 != 0)).any().any():
    p.error("states must be integer -1, 0 or positive")
tree = cas.data.CassiopeiaTree(character_matrix=matrix)
cas.solver.VanillaGreedySolver().solve(tree, collapse_mutationless_edges=True)
Path(args.output_newick).write_text(tree.get_newick().rstrip() + "\n", encoding="utf-8")
