# Popcorn: Trans-Ancestry rg
**Goal:** Quantify shared causal architecture between two ancestries (e.g. EUR T2D vs EAS T2D) under a MAF-LD-aware model.

**Approach:** Compute cross-population LD scores once per ancestry pair; fit rg using sumstats from each population.

```bash
# Step 1: cross-population LD scores (one-time per ancestry pair)
popcorn compute -v 1 \
    --bfile1 1kg_EUR \
    --bfile2 1kg_EAS \
    --SNPs_to_store 20000 \
    --gen_effect \
    eur_eas_scores.txt

# Step 2: fit cross-population rg
popcorn fit -v 1 \
    --cfile eur_eas_scores.txt \
    --gen_effect \
    --sfile1 t2d_eur.sumstats.txt \
    --sfile2 t2d_eas.sumstats.txt \
    t2d_eur_eas_rg.txt
```

`--gen_effect` (which must be passed at BOTH `compute` and `fit`) selects the genetic-effect model and reports `rho_ge` (correlation of causal effect sizes); omitting it at both steps yields `rho_gi` (correlation of variant-level impacts, MAF-weighted). A single run reports one or the other, so run both modes to report both. When MAFs differ markedly across populations, `rho_ge` and `rho_gi` diverge; both are biologically meaningful and report-worthy. Effective N per population must be > ~5000 for stable estimates.
