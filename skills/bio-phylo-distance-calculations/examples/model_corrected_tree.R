# Model-corrected distance tree in ape. Deliverable: a TN93+gamma distance matrix and a FastME
# balanced-minimum-evolution tree with bootstrap -- the correction is where the biology lives,
# unlike Biopython's identity-only DistanceCalculator. Outputs go to tempdir(), no strays.
# Runs on data/sim8.fasta (8 taxa, 300bp, simulated under a known tree, data/sim8_true.nwk) and
# reports the Robinson-Foulds distance to that true tree, so this checks against a known answer
# instead of a handful of toy strings -- and contrasts with build_nj_tree.py on the same data.
# Reference: ape 5.8+ | Verify API if version differs

library(ape)
library(phangorn)

# Run from the examples/ directory (or adjust the path): loads the shipped simulated alignment.
aln <- read.dna('data/sim8.fasta', format = 'fasta')

# TN93 = two transition rates + unequal base freqs; gamma applies ASRV (alpha < 1 = strong heterogeneity).
# Estimate alpha by ML (phangorn) instead of assuming a value -- same as the SKILL.md snippet.
aln_phy <- phyDat(aln, type = 'DNA')
fit <- optim.pml(pml(nj(dist.dna(aln, model = 'JC69')), aln_phy, k = 4), optGamma = TRUE,
                  model = 'GTR', rearrangement = 'none')
alpha <- fit$shape

# FastME balanced minimum evolution: the modern best distance tree (searches, not a single greedy pass).
# Define the method ONCE so the bootstrap replicates use exactly the correction that built the tree.
method <- function(x) fastme.bal(dist.dna(x, model = 'TN93', gamma = alpha), nni = TRUE, spr = TRUE)
tree <- method(aln)

# Bootstrap: 100-1000 reps standard; this is sampling PRECISION, not accuracy.
set.seed(42)
bs <- boot.phylo(tree, as.matrix(aln), method, B = 100, quiet = TRUE)
tree$node.label <- bs                               # attach support (replicate counts out of B)
cat('Bootstrap counts per internal node (of 100):', bs, '\n')

# Check against the known-true topology (dist.topo = Robinson-Foulds split distance):
true_tree <- read.tree('data/sim8_true.nwk')
rf <- dist.topo(unroot(tree), unroot(true_tree))
cat('Robinson-Foulds distance to the true tree:', rf, '(0 = identical unrooted topology)\n')
cat('(contrast with build_nj_tree.py on the same data, which uses uncorrected identity distance)\n')

out <- file.path(tempdir(), 'distance_tree.nwk')
write.tree(tree, out)
cat('Wrote model-corrected FastME tree to', out, '\n')
print(tree)
