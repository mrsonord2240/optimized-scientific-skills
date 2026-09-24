# Genome-Wide Library Catalog and PAM Variants

## Genome-Wide Library Selection

| Library | Year | Modality | Size (genes x guides) | sgRNA rules | Notable |
|---------|------|----------|-----------------------|-------------|---------|
| GeCKOv2 | 2014 | Cas9 KO | ~19k x 6 (~123k) | Exon position + off-target specificity (predates Rule Set 1) | Older; legacy datasets still use it |
| Avana | 2016 | Cas9 KO | 110,257 as published; DepMap screens a ~4-guide subset (Meyers 2017: 70,086 after filtering, 17,670 genes) | Rule Set 1 | Still the Broad's primary Cas9 library; CERES->Chronos changed in 2021, not the library |
| Brunello | 2016 | Cas9 KO | ~19k x 4 (~77k) | Rule Set 2 + CFD | Modern standard for new screens |
| TKOv3 | 2017 | Cas9 KO | ~18k x 4 (~71k) | Hart on/off-target | Bagel/BAGEL2-optimized |
| Humagne | 2020 | enAsCas12a | ~19.8k x 1 dual-guide construct (~20k per set) | enAsCas12a rules | Compact Cas12a sets C and D |
| Horlbeck CRISPRi v2 | 2016 | dCas9-KRAB | ~18k x 5 (~104k) | Horlbeck CRISPRi rules | First-gen, still widely used |
| Dolcetto | 2018 | dCas9-KRAB | ~19k x 3 per set (114,061 across Sets A+B) | Horlbeck + Rule Set 2 | Modern CRISPRi standard |
| Horlbeck CRISPRa | 2016 | dCas9-VP64 | ~18k x 5 (~104k) | Horlbeck CRISPRa rules | Original CRISPRa |
| Calabrese | 2018 | dCas9-VP64 | ~18.9k x 3 per set (113,238 across Sets A+B) | Tight TSS window | Modern CRISPRa standard |
| Inzolia | 2024 | enAsCas12a | ~49k arrays: 19,687 genes (2 arrays each) plus ~4,435 paralog pairs | enAsCas12a rules | Paralog-pair multiplex; ~30% smaller than a typical Cas9 library |
| in4mer | 2024 | Cas12a (4-guide) | Custom | enAsCas12a multiplex | Triple/quadruple KO per cassette |

**dAUC trajectory (essentiality benchmark):** GeCKOv2 < Avana < Brunello/TKOv3 (Doench 2016 + Hart 2017). Moving from 4 to 6 sgRNAs/gene gives diminishing returns; the larger gain is moving from Rule Set 1 to Rule Set 2.

**Cost-coverage tradeoff:** A 77k-guide Brunello at 500x cells/sgRNA needs 38.5M cells in pool, scalable. A 117k-guide Calabrese at 500x needs 59M cells -- often the deciding factor against CRISPRa for difficult-to-grow lines.

## PAM Variants and Alternative Cas Enzymes

| Enzyme | PAM | Spacer length | Best for |
|--------|-----|---------------|----------|
| SpCas9 (WT) | NGG | 20 nt | Standard pooled screens; broadest library support |
| eSpCas9, SpCas9-HF1 | NGG | 20 nt | Lower off-target rate; use for therapeutic-grade nomination |
| SpCas9-NG | NG | 20 nt | Expanded targeting (~4x coverage); accept lower activity per guide |
| SpRY | NRN / NYN | 20 nt | Near-PAMless; coverage at every position; ~50% lower per-guide activity |
| SaCas9 | NNGRRT | 21 nt | AAV-packageable (small ORF); rarely used in pooled screens |
| AsCas12a, LbCas12a | TTTV | 23 nt | AT-rich regions; staggered cut; lower expression noise |
| enAsCas12a (DeWeirdt 2021) | Expanded TTTV + several non-canonical | 23 nt | Combinatorial / paralog screens |

**Decision rule:** If the screen requires every possible TSS position (saturation tiling, dense regulatory dissection), use SpRY despite lower activity; otherwise, NGG is best because the on-target predictors were trained on it.
