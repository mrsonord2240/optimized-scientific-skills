# Capra-Singh Jensen-Shannon divergence

Read when scoring catalytic or functional residues on a protein MSA with JSD. Assumes the normalised `alignment` and the imports from "Required Import and Normalisation" in `SKILL.md`.

### Capra-Singh Jensen-Shannon Divergence

**Goal:** Score columns by divergence from a residue-frequency background, with a window-smoothed neighbour penalty for catalytic residue prediction.

**Approach:** Compute JSD between the column distribution (letters outside the background ignored) and a residue-frequency background (Capra & Singh 2007 Bioinf used BLOSUM62-derived; the example uses Robinson & Robinson 1991 PNAS, which gives effectively-equivalent column ranking: Spearman 0.98 and top-10 overlap 9/10 against the authors' script, which also applies Henikoff sequence weights that this example omits), then mix with the windowed-neighbour mean. Defaults `window=3`, `lambda_window=0.5` track catalytic-residue annotation in the Catalytic Site Atlas. Full implementation: `examples/capra_singh_jsd.py`.

```python
def capra_singh_score(alignment, background=None, window=3, lambda_window=0.5):
    # raw[i] = JSD(column_i, background) * (1 - gap_fraction_i)   # gap-penalty per Capra-Singh reference impl
    # smoothed[i] = (1 - lambda) * raw[i] + lambda * mean(raw[i-window:i+window+1] excluding i)
    ...
```

**Threshold for catalytic-residue prediction:** Capra & Singh 2007 (Bioinf 23:1875) report AUC ~0.94 and Top-30 score ~0.75 on the Catalytic Site Atlas using JSD with neighbor mixing (window=3, lambda=0.5); the paper does not prescribe a single threshold. Choose by ROC tradeoff for the specific use case. ConSurf-derived rate4site (rate-of-evolution) is competitive but requires a phylogenetic tree; Capra-Singh JSD is the alignment-only equivalent.
