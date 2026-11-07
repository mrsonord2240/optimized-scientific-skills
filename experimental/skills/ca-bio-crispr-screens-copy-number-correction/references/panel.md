## Chronos (Dempster 2021) - Joint Population-Dynamics + CN Model

**Goal:** Estimate gene fitness while jointly accounting for copy-number-driven depletion, screen quality, and longitudinal cell-population dynamics.

**Approach:** Model the cell population over time as an ODE driven by per-gene fitness effects; add a separate term for copy-number-driven depletion; estimate all parameters via maximum-likelihood with regularization. Outputs a "gene effect score" normalized against the empirical distributions of essential and non-essential reference genes.

```python
# Chronos (pip install crispr_chronos, or pip install git+https://github.com/broadinstitute/chronos)
import numpy as np
import chronos
from chronos.hit_calling import get_probability_dependent

# Chronos initializes model parameters with NumPy draws. Seed *before* constructing the model,
# and record the Chronos version, seed, input hashes, epochs, and platform in the run manifest.
# This makes a CPU rerun materially comparable; do not promise bitwise-identical estimates across
# TensorFlow/BLAS versions or hardware. Re-run important panels and report their stability.
np.random.seed(20260924)  # replace with and record a project-specific seed

# Every input is a dict keyed by library name, not a bare DataFrame, and Chronos checks the
# schema below (chronos.check_inputs) before training:
#
# 1. readcounts    rows = sequence_ID (one per sequenced sample), columns = sgRNA.
#                  This is the transpose of a MAGeCK count table; Chronos rejects the other
#                  orientation with "Chronos expects readcounts to have guides as columns,
#                  sequence IDs as rows. Is your data transposed?"
# 2. sequence_map  columns sequence_ID, cell_line_name, days, pDNA_batch. Plasmid samples get
#                  cell_line_name 'pDNA' and days 0, and every pDNA_batch used by a late
#                  timepoint needs a pDNA row in the same library.
# 3. guide_gene_map  columns sgrna, gene; one gene per sgRNA (duplicated sgRNAs are rejected).
# 4. negative_control_sgrnas  dict library -> sgRNAs of non-targeting or known non-essential
#                  genes. Optional in the signature only: without it, chronos.Chronos(...)
#                  construction itself raises ValueError('excess_variance was passed as dict
#                  without key for <library>: {}') -- before train() is ever called. An
#                  UnboundLocalError for 'prior_variance' does occur one frame deeper, inside
#                  _estimate_excess_variance, but Chronos catches it internally and re-raises
#                  the ValueError above every time; it never reaches the caller under default
#                  arguments (checked on Chronos 2.3.15).
# 5. Copy-number profile: cell lines x genes, applied AFTER training (see below).

chronos.check_inputs(                      # fail fast on schema and orientation
    readcounts={'screen': readcounts_df},  # rows = sequence_ID, columns = sgRNA
    guide_gene_map={'screen': guide_gene_map},
    sequence_map={'screen': sequence_map},
)

model = chronos.Chronos(
    sequence_map={'screen': sequence_map},
    guide_gene_map={'screen': guide_gene_map},
    readcounts={'screen': readcounts_df},
    negative_control_sgrnas={'screen': negative_control_sgrnas},
)
model.train(nepochs=301)
gene_effects = model.gene_effect                      # attribute, not a method call; lines x genes

# Copy-number correction is a separate post-hoc step, not a constructor argument. It needs at
# least 3 cell lines (RuntimeError below that), the CN frame must be lines x genes and cover
# every gene in gene_effects, and it returns two objects: the corrected matrix and the
# per-gene CN shifts it fitted.
gene_effects_cn, cn_shifts = chronos.alternate_CN(gene_effects, copy_number_df)
gene_probabilities = get_probability_dependent(gene_effects_cn, negative_control_genes, positive_control_genes)
```

**DepMap convention:** A gene-effect score <-1 corresponds to "essential" in that cell line; <-0.5 is "depleting." Each DepMap release (quarterly) provides Chronos gene effects and probabilities.

**Critical:** Chronos benefits most from longitudinal data (multiple timepoints per cell line) but can run with multiple cell lines at a single timepoint, and it trains on a single line with replicates too. Copy number is optional to *training*: Chronos trains without it and `alternate_CN` applies the correction afterwards. The correction itself has a hard minimum of 3 cell lines, so for one or two lines, CN-correct with CRISPRcleanR (pre-hoc, on the counts) and use Chronos only for the gene-effect scores, or not at all.

## CERES (Legacy, Superseded by Chronos)

**Meyers RM et al 2017 *Nat Genet* 49:1779** introduced the first formal CN-correction method at DepMap scale. CERES decomposes per-sgRNA LFC as `sgRNA_efficacy * gene_effect - CN_term(copy_number)`, fitting jointly. Superseded by Chronos at DepMap in 2021 due to Chronos' better handling of screen quality and longitudinal data. CERES remains useful for cross-validation.
