# Background

Not needed to run a route.

A ctDNA assay is a molecule-counting experiment at the Poisson edge. The mutant signal is a fixed, tiny number of physical templates, and the limit of detection is set by two ceilings: how many genome equivalents were sampled, and how low the background error floor was driven. Under the 3.3 pg per haploid genome convention, 1 ng holds 303.03 GE. Past sampling saturation (a dedup plateau in UMI families) extra read depth re-sequences the same molecules and adds no information.

Definitions (CLSI EP17-A2): LoB = mean_blank + 1.645 SD_blank; LoD = LoB + 1.645 SD_low (about 95% detection at LoD); LoD95 is the point on a probit or logistic detection curve where detection probability is 95%; LoQ is the lowest level measurable with stated precision and is at least LoD. Per-locus LoD is typically 0.05-0.1% VAF.

Binomial multi-locus integration over equal-VAF, independent Poisson presence probabilities is a lower bound, not an assay LoD95. The fit script uses the Gaussian LoB shortcut and a probit model on log10 VAF.
