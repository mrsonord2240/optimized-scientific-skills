import argparse
import numpy as np
from Bio import AlignIO, Phylo
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor

parser = argparse.ArgumentParser()
parser.add_argument('alignment')
parser.add_argument('output')
args = parser.parse_args()
aln = AlignIO.read(args.alignment, 'fasta')
arr = np.array([list(rec.seq) for rec in aln])
keep = (arr != '-').all(axis=0)
if not keep.any():
    parser.error('no alignment columns remain after gap filtering')
aln = MultipleSeqAlignment([SeqRecord(Seq(''.join(row[keep])), id=rec.id) for row, rec in zip(arr, aln)])
matrix = DistanceCalculator('identity').get_distance(aln)
tree = DistanceTreeConstructor().nj(matrix)
Phylo.write(tree, args.output, 'newick')
print(matrix)
