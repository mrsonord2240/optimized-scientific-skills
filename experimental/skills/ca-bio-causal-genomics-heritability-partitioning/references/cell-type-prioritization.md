# Cell-type prioritization, Finucane 2018 (reference for heritability-partitioning SKILL.md)

## Cell-Type Prioritization (Finucane 2018)

Finucane 2018 Nat Genet 50:621 introduced cell-type-specific S-LDSC: partition heritability against ENCODE / Roadmap chromatin marks (H3K4me3, H3K27ac, H3K4me1, DNase, ATAC) tissue-by-tissue, retain per-tissue p-value adjusting for the baseline model. Trait-relevant tissue = top-ranked tissue with p < 0.05/N_tissues (Bonferroni for ~200 tissues, threshold ~2.5e-4).

**Goal:** Rank tissues / cell types by their per-annotation contribution to trait heritability.

**Approach:** Build per-cell-type LD scores from chromatin-marker BED files; compile `.ldcts` manifest (one row per cell type: name, ldscore prefix, control ldscore); run `--h2-cts` and interpret per-cell-type coefficient p-value.

```bash
# Cell-type prioritization example workflow (Finucane 2018)
# Inputs: trait.sumstats.gz (munged), <cts>.ldcts manifest, baseline annotations,
#         eur_w_ld weights, 1000G EUR frequency files

ldsc.py \
    --h2-cts trait.sumstats.gz \
    --ref-ld-chr 1000G_EUR_Phase3_baseline/baseline. \
    --ref-ld-chr-cts Multi_tissue_chromatin.ldcts \
    --w-ld-chr weights_hm3_no_hla/weights. \
    --out trait_cts
# trait_cts.cell_type_results.txt: Name, Coefficient, Coefficient_std_error, Coefficient_P_value
# Apply Bonferroni at 0.05 / nrow; top tissues are trait-relevant
```

The `.ldcts` manifest is tab-separated, one row per cell type: `<name>\t<ldscore_prefix>,<control_ldscore_prefix>` (the control prefix is optional -- omit the comma if there is none). `--ref-ld-chr-cts` and the baseline `--ref-ld-chr` both need chromosome-split LD score files (e.g. `prefix1.l2.ldscore.gz` ... `prefix22.l2.ldscore.gz`), not a single-file prefix.

Published `.ldcts` files cover GTEx tissues, Roadmap epigenome, immune cell types, and scATAC clusters. Custom .ldcts for novel cell types requires computing per-cell-type LD scores from a chromatin BED via `ldsc.py --l2 --bfile ... --annot <cell>.annot.gz`.
