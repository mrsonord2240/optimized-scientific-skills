sim8.fasta / sim8_true.nwk: synthetic 8-taxon, 300bp DNA alignment simulated with
IQ-TREE 2.4.0 AliSim (--alisim, -m JC69, --seqtype DNA, --length 300, --seed 42) along
sim8_true.nwk (branch lengths 0.04-0.90, a mix of short and long branches so distance
correction and saturation both matter). Used by build_nj_tree.py and
model_corrected_tree.R so the examples recover a tree against a known answer instead of
a 12-20bp toy string.
