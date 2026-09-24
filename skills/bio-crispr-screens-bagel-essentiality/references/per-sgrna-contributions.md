# Per-sgRNA Bayes Factor contributions

> Moved out of SKILL.md. "Below" and "Failure Modes" in the text refer to sections of SKILL.md.

## Bayesian Reasoning Per Sgrna

**Why this matters:** BAGEL2 computes per-sgRNA contributions; a gene with 4 sgRNAs each contributing +5 to BF gets +20 total. A gene with 3 sgRNAs contributing +5 and 1 sgRNA contributing -3 (off-target or low-efficacy) gets +12 net.

**Approach:** Add `-r/--sgrna-bayes-factors` to the `bf` command (found via `BAGEL.py bf --help`;
required for per-sgRNA output).

```bash
BAGEL.py bf \
    -i foldchange.foldchange \
    -o bayes_factor_sgrna.txt \
    -e CEGv2.txt -n NEGv1.txt \
    -c Sample1,Sample2,Sample3 \
    -s 42 -r                               # -r: per-sgRNA BF contributions
# Output columns: RNA  GENE  <sample columns>  BF  (one row per sgRNA)
```

Verified on real HAP1 TKOv3 data: for RPS3 (the example gene below), the 4 per-sgRNA
BF values summed to 78.9 against a gene-level BF of 80.8 (within 2.4%), confirming the
additive-LLR model this section describes.

**Critical:** When per-sgRNA contributions are very heterogeneous (one sgRNA dominates BF), the gene is "guide-of-one"; verify with JACKS efficiency analysis or apply the second-best-sgRNA rule from [[hit-calling]].
