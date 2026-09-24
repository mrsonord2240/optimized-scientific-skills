# JACKS efficacy prior and guide diagnostics

## Build Library-Wide Efficacy Prior from Reference Screens

**Goal:** Transfer learned efficacy from a large public screen panel to a new small screen.

**Approach:** Run JACKS on the reference panel (e.g. DepMap CRISPR screens with TKOv3 or Brunello), extract per-sgRNA efficacy posterior, and supply it as the prior for a new screen.

```python
def extract_efficacy_prior(reference_jacks_results):
    '''Build per-sgRNA efficacy prior (mean + std) from a large reference screen.'''
    df = pd.read_csv(reference_jacks_results, sep='\t')
    prior = df[['sgrna', 'X1', 'X2']]      # --reffile requires these exact column names; do not rename
    return prior

# Use in new JACKS run via --reffile <path>
# Reference: Allen 2019 Genome Research 29:464; efficacy-aware testing enables ~2.5x smaller screens (fewer replicates/guides),
# but only when the reference is the same library and a similar cell context. Reference panels are ~50 cell lines, ~10k screen days.
```

## Per-sgRNA Efficacy Diagnostics

**Goal:** Identify low-efficacy guides for library refinement.

**Approach:** Examine the distribution of inferred efficacies; guides below 0.3 are likely non-functional and should be excluded from re-designed libraries.

```bash
python scripts/efficacy_summary.py jacks_out_grna_JACKS_results.txt guidemap.txt --low 0.3 --out low_eff_by_gene.tsv
# or: from efficacy_summary import efficacy_summary; summary, by_gene = efficacy_summary(grna_path, guidemap_path)
```

The grna file has only `sgrna`, `X1`, `X2`, so genes come from the guide map; the script raises `ValueError` if any result guide is absent from the map (naming mismatch). It reports total guides, low-efficacy count and percentage, median and quartiles of `X1`, and the per-gene fraction of low-efficacy guides.

**Critical:** Genes where every guide is low-efficacy will show no signal regardless of biology. Filter from interpretation; flag for re-design with updated rules (Brunello / TKOv3). For a v2 library, drop the bottom 25% of guides by efficacy; every gene should end with all guides at efficacy >0.4 (Brunello v2 / Avana v2 convention).
