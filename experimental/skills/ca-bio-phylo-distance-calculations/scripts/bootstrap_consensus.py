import argparse
import random
from Bio import AlignIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Phylo.Consensus import bootstrap_trees, majority_consensus

parser = argparse.ArgumentParser()
parser.add_argument('alignment')
parser.add_argument('outgroup')
parser.add_argument('output')
parser.add_argument('--replicates', type=int, default=50)
args = parser.parse_args()
alignment = AlignIO.read(args.alignment, 'fasta')
random.seed(42)
constructor = DistanceTreeConstructor(DistanceCalculator('identity'), 'nj')
replicates = []
for tree in bootstrap_trees(alignment, args.replicates, constructor):
    tree.root_with_outgroup({'name': args.outgroup})
    replicates.append(tree)
Phylo.write(majority_consensus(replicates, cutoff=0.5), args.output, 'newick')
