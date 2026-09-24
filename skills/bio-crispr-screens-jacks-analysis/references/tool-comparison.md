# JACKS vs MAGeCK and BAGEL2

## Comparing JACKS, MAGeCK, BAGEL2

| Property | JACKS | MAGeCK | BAGEL2 |
|----------|-------|--------|--------|
| Statistical framework | Variational Bayes | NB GLM + alpha-RRA / MLE | Bayes factor on per-sgRNA fold change |
| Models guide efficacy | Yes (jointly) | No (optional fixed input) | No |
| Multi-screen joint | Yes (native) | Limited (MLE design matrix) | No (per-screen) |
| Speed | Slow (variational inference) | Fast | Fast |
| Output | gene effect + sgRNA efficacy | beta or RRA score | Bayes Factor |
| Best for | Multi-screen joint analyses, library calibration | General-purpose, single screen | Essentiality classification |

**Reconciliation:** Hits identified by JACKS AND MAGeCK are high confidence. JACKS-only hits typically reflect strong gene signals where one or two guides were dragging down MAGeCK; verify the up-weighted high-efficacy guides have the expected sign. MAGeCK-only hits at FDR <0.05 may be single-guide outliers; check sgrna_summary for guide-level dispersion.

## Reconciliation: When JACKS and Other Tools Disagree

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| JACKS significant, MAGeCK not | One low-efficacy guide dragged MAGeCK; JACKS down-weighted it | Trust JACKS if 3+ high-efficacy guides agree |
| MAGeCK significant, JACKS not | All guides have similar efficacy; JACKS prior shrinks signal | Verify per-guide LFC consistency in MAGeCK sgrna_summary |
| JACKS efficacy ~0.5 for all guides | Efficacy prior mismatched to the chemistry, or weak signal overall (`--apply_w_hp` acts on gene effects, not efficacy) | Confirm chemistry and library match; see "Efficacy collapsed near zero" in `failure-modes.md` for a matched `--reffile` or prior override |
| Gene effect different sign from MAGeCK | Multi-screen pooling created mean effect different from single-screen | Run per-screen separately to confirm |
