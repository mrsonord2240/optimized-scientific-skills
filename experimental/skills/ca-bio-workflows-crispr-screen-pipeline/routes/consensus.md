# Tier consensus across methods

```bash
python scripts/consensus.py tier_consensus.csv mageck=essentiality_rra.gene_summary.txt \
    bagel=bayes_factor.txt drugz=drugz_output.txt
```

Give at least two of `mageck=`, `bagel=`, `drugz=`, all from depletion results on the same count table and baseline. Hit rules: MAGeCK `neg|fdr < 0.05`, BAGEL2 `BF > 6`, drugZ `fdr_synth < 0.05` (`fdr=` and `bf=` override).

- Tier-1: all three agree; Tier-2: two; Tier-3: one, exploratory. With two methods the best tier is Tier-2.
- The BAGEL2 file must come from a seeded run (`routes/bagel2.md`).
- For a high-stakes target nomination, report Tier-1 and Tier-2 only.
- Validate Tier-1 orthogonally.

Done when `tier_consensus.csv` exists and tier counts are reported.
