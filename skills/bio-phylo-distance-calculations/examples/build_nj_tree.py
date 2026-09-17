'''Build a Neighbor Joining tree from an alignment. Deliverable: a fast NJ tree -- but the
distance here is IDENTITY-ONLY (a p-distance), NOT a model correction. For divergent DNA use
ape dist.dna(model='TN93') for a model-corrected matrix; the matrix is the limitation, not NJ.
Runs on data/sim8.fasta (8 taxa, 300bp, simulated under a known tree) and checks the recovered
topology against data/sim8_true.nwk with the Robinson-Foulds distance, so the example tests
against a known answer instead of a handful of toy strings.'''
# Reference: biopython 1.83+, scikit-bio 0.7+ | Verify API if version differs

import io
import os

from Bio import AlignIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from skbio import TreeNode

data_dir = os.path.join(os.path.dirname(__file__), 'data')
alignment = AlignIO.read(os.path.join(data_dir, 'sim8.fasta'), 'fasta')

calculator = DistanceCalculator('identity')   # identity-only: p-distance, NO multiple-hit correction
dm = calculator.get_distance(alignment)        # divergent data needs ape dist.dna for JC/K80/TN93/LogDet

print('Distance Matrix (uncorrected p-distance):')
print(dm)

constructor = DistanceTreeConstructor(calculator, 'nj')
tree = constructor.build_tree(alignment)       # NJ is consistent ONLY if the distances are correct/additive
tree.ladderize()

print('\nNeighbor Joining Tree:')
Phylo.draw_ascii(tree)

# Check against the known-true topology (Robinson-Foulds, as a proportion of the max possible):
handle = io.StringIO()
Phylo.write(tree, handle, 'newick')
estimated = TreeNode.read(io.StringIO(handle.getvalue()))
true_tree = TreeNode.read(os.path.join(data_dir, 'sim8_true.nwk'))
rf = estimated.compare_rfd(true_tree, proportion=True)
print(f'\nRobinson-Foulds distance to the true tree (proportion): {rf:.3f}')
print('(identity distance is uncorrected, so some topological error here is expected;')
print(' contrast with model_corrected_tree.R on the same data.)')
