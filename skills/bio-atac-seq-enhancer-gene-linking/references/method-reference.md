# Enhancer-gene linking method reference

Executed here (2026-09-30, `scripts/run_abc.sh` and the pipelines on the ABC chr22 K562 example): ABC-Enhancer-Gene-Prediction v1.1.2 and `main` 92ac503 (`workflow/scripts/` layout; ABC has no PyPI version, only git tags, and `src/` is not used in any tag checked); ENCODE-rE2G `main` d039062 (`v1.0.0`, `v2.0.0igvf` are its tags); python 3.10, pandas 2.3.3, pyranges 0.1.4, hic-straw 1.3.1, MACS2 2.2.9.1 (needs a glibc-compat shim for `__log_finite` on glibc >= 2.31) or MACS3 3.0.4, bedtools 2.31.1 (ENCODE-rE2G pins 2.29.2 and scikit-learn 1.2.1), snakemake 7 with `pulp<2.8`. Origin-stated floors not executed: Cicero 1.20+, GenomicInteractions 1.36+, FitHiChIP 9.1+, HiC-Pro 3.1+, FAN-C 0.9+. Verify with `<tool> --version`, `packageVersion('<pkg>')`, or `pip show <package>`.

## Algorithmic taxonomy

| Method | Inputs | Mathematics | Strength | Fails when |
|--------|--------|-------------|----------|------------|
| ABC (Fulco 2019, Nasser 2021) | ATAC + H3K27ac + Hi-C/Micro-C | ABC = (Activity_E x Contact_E,G) / sum_e(Activity_e x Contact_e,G); threshold calibrated per input (0.012-0.027) | Mechanistically grounded; published gold-standard for human cell lines | Best with matched Hi-C / Micro-C; without it, contact is a distance powerlaw or a cell-type-average Hi-C, both lower performing per the ABC docs |
| ENCODE-rE2G (Gschwind 2026; preprint 2023) | DNase or ATAC, optional H3K27ac, intact or megamap Hi-C | Logistic regression on ABC-derived features plus distance and promoter/enhancer context, trained on K562 CRISPR-validated links | ENCODE 4 standard; thresholds set for 70% recall on the CRISPR benchmark | No powerlaw or average-Hi-C model; models are selected by input type, not cell type; training needs matching CRISPR data |
| Cicero (Pliner 2018) | scATAC peak-cell matrix | Graphical lasso on metacell co-accessibility | ATAC-only; works without Hi-C | Less concordant with Hi-C than ABC; cis-distance-limited; alpha-sensitive |
| HiChIP H3K27ac + FitHiChIP | H3K27ac HiChIP | Statistically significant loops at FDR < 0.05 | Direct experimental loop measurement; cell-type-specific; orthogonal to ATAC | Requires HiChIP wet-lab; only captures loops within HiChIP resolution (~10 kb) |
| Hi-C + HiCCUPS | Bulk Hi-C | Fold-enrichment loop calling | Most-validated 3D contact method | Resolution typically 5-25 kb; misses sub-loop fine structure |
| Capture Hi-C / PCHi-C (CHiCAGO) | Promoter Capture Hi-C | Asymptotic CHiCAGO score | High-resolution promoter-anchored | Wet-lab cost; promoter capture only |
| EpiMap (Boix 2021) reference | None (pre-computed lookup) | Bulk-derived enhancer-gene predictions in 833 epigenomes | Fast, comprehensive | Cell-type-agnostic for tissues outside the reference set |
| GeneHancer / FANTOM5 (legacy) | None (pre-computed lookup) | Pre-computed; varied methods per database | Comprehensive lookup; widely cited | Older; less reliable than ABC for cell-type-specific |

Methodology evolves; verify against current Engreitz lab releases (ABC), ENCODE 4 publications (ENCODE-rE2G), and Mumbach 2017 (HiChIP) before locking pipelines.

## ABC mathematics

For each candidate (enhancer E, gene G) pair within the cis window (default 5 Mb):

```
ABC(E -> G) = Activity_E * Contact_E,G / sum_{all e in window}(Activity_e * Contact_e,G)
```

- **Activity_E** = ATAC reads at E * H3K27ac reads at E (geometric mean of normalized signals; reflects "enhancer strength")
- **Contact_E,G** = Hi-C/Micro-C contact frequency from E to G's TSS (after distance-correction)
- **Window** = +/- 5 Mb cis (`predict.py --window`, default 5000000)

The threshold is calibrated per input combination in ABC's `reference/abc_thresholds.tsv` (accessibility type, H3K27ac presence, contact type: 0.012-0.027); 0.02 is only the fallback the pipeline uses when no row matches. The Snakemake pipeline maps `avg` Hi-C to `avg_hic`, which has no row, so it falls back to 0.02; `run_abc.sh` uses the table's `avg` rows. Lower cut-offs than the calibrated value are exploratory.

Contact options, best to worst per the ABC docs (highest performance with cell-type-specific Hi-C): matched Hi-C; a powerlaw of genomic distance (Fulco 2019; explains >70% of Hi-C variance, and the only choice for non-human organisms); or a cell-type-average Hi-C map (Nasser 2021), which adds cell-type-invariant features such as CTCF loops but needs the 58 GB ENCODE ENCFF134PUN bed, extracted per chromosome with `extract_avg_hic.py` (1-2 h). Either fallback is a documented degradation to report.

## ABC standard pipeline notes

`scripts/run_abc.sh` implements the steps. Points to know:

- `run.neighborhoods.py` counts reads directly from the accessibility and H3K27ac BAMs; no bigWig is used, so the script makes none.
- Candidate regions come from `makeCandidateRegions.py` on MACS2 summit peaks (`--call-summits` is required). The TSS regions are included on purpose: elements within 500 bp of a TSS are classed `promoter`, and `filter_predictions.py` drops them unless they are the gene's own promoter. Removing promoters from the candidates beforehand would drop the promoters' Activity x Contact from each gene's normalising sum and the score of 1 given to a gene's own promoter.
- `predict.py` requires `--chrom_sizes`, `--accessibility_feature` (`DHS` or `ATAC`), `--hic_gamma`/`--hic_scale` (positive powerlaw parameters; ABC's config carries the average-Hi-C values), and `--hic_pseudocount_distance` (required on `main`, default 1e6 in v1.1.2). `--hic_type` is `hic` (`.hic` file or URL) or `juicebox | bedpe | avg` (directory); the Snakemake pipeline scores the no-Hi-C case with `--score_column powerlaw.Score`.
- `predict.py` writes `EnhancerPredictionsAllPutative.tsv.gz` (expressed genes) and `...AllPutativeNonExpressedGenes.tsv.gz`; `filter_predictions.py` applies the threshold, removes non-self promoters, and writes the full, slim, bedpe, and gene-stats tables.

## ENCODE-rE2G differences from ABC

ENCODE-rE2G (Gschwind et al 2026; the README also links the 2023 preprint) is a logistic regression built on ABC:

- **Training:** CRISPR-validated links; the repository states training is for K562-matched CRISPR data only.
- **Features:** per-model `feature_table.tsv` (ABC score and contact, distance, promoter activity, nearby-enhancer activity, ubiquitous-gene flag, candidate counts).
- **Model selection is automatic** from the biosample row: `{dhs|atac}[_h3k27ac]_{intact_hic|megamap}` in `models/`, plus an `extended` model chosen only via a `model_dir` column. `HiC_type=avg` or an empty `HiC_file` finds no model directory and the pipeline raises an exception. A cell type without its own Hi-C uses the megamap `.hic` (`MEGAMAP_HIC_FILE` in `config/config.yaml`).
- **Output:** `ENCODE-rE2G.Score` in [0,1] per pair; `encode_e2g_predictions_threshold{t}.tsv.gz` applies the model's `threshold_*` file (0.179-0.298).

Use ENCODE-rE2G when cell-type or megamap Hi-C is available; use ABC powerlaw when no Hi-C is wanted. Because models are keyed by input type, cell-type proxies apply to Hi-C only: name the Hi-C source. Compare the two methods on the same candidate regions.

## Per-tool failure modes

### ABC: wrong cell-type-matched Hi-C

**Trigger:** Using K562 Hi-C contact when actual cell type is GM12878.
**Mechanism:** Contact frequencies differ across cell types at compartment and TAD boundaries.
**Symptom:** ABC predictions concentrate at known K562-specific loci even when ATAC data is from GM12878.
**Fix:** Use cell-type-matched Hi-C or Micro-C. If unavailable, use the powerlaw or ABC's cell-type-average Hi-C, both documented fallbacks with acknowledged degradation. Document the proxy in methods. Hi-C resolution matters too: ABC uses 5 kb and states that other resolutions change model behaviour.

### ABC: H3K27ac normalization

**Trigger:** H3K27ac ChIP-seq with different sequencing depth than ATAC.
**Mechanism:** Activity is the geometric mean of accessibility and H3K27ac read counts over the candidate; without normalization sequencing depth and assay differences shift the score distribution away from the one the thresholds were calibrated on.
**Symptom:** Activity skewed to one assay, or scores that do not match the calibrated thresholds.
**Fix:** Keep `--qnorm reference/EnhancersQNormRef.K562.txt` (ABC ships only the K562 reference; `run_abc.sh` applies it by default, `QNORM=none` disables it, and ABC's config warns that disabling it gives scores that do not follow the preset thresholds).

### ABC: promoter links reported as enhancers

**Symptom:** Predictions concentrate at TSSs.
**Fix:** Keep the TSS regions in the candidate set and drop `promoter`-class links after scoring (`filter_predictions.py`); see the pipeline notes above.

### ENCODE-rE2G: input combination without a model

**Trigger:** A biosample with no Hi-C (powerlaw) or `HiC_type=avg`.
**Mechanism:** No such model is shipped, so model selection raises an exception.
**Fix:** Provide cell-type Hi-C or point `HiC_file` at the megamap `.hic` (a documented Hi-C proxy to report). The coefficients were fit on K562 CRISPR data, so treat other cell types as extrapolation. Custom training needs matching CRISPR data.

### Cicero: no Hi-C concordance benchmark

**Trigger:** Reporting Cicero connections as enhancer-gene calls without external validation.
**Mechanism:** Cicero is statistical co-accessibility; correlation with Hi-C contacts is ~30-50%. Many strong Cicero connections are not Hi-C-validated.
**Fix:** When Hi-C is available, cross-validate and report both (switch to ABC). When only ATAC, use Cicero with the explicit caveat that connections are co-accessibility hypotheses, not contact predictions.

### HiChIP: loop-calling threshold

**Trigger:** Default FitHiChIP at FDR < 0.05.
**Mechanism:** HiChIP loops are abundant (10k-100k per dataset); FDR alone leaves a long tail of weak loops.
**Fix:** Threshold at FDR < 0.05 and contacts per loop >= 5; or use the top N most significant where N is the expected loop count for the cell type.

### EpiMap / GeneHancer / FANTOM5: cell-type-agnostic limitation

**Mechanism:** These references aggregate across many tissues and experiments; cell-type-specific connections are diluted.
**Fix:** Use as baseline or sanity check only. ABC or ENCODE-rE2G in the actual cell type is preferred.

## CRISPRi-FlowFISH validation

CRISPRi-FlowFISH (Fulco 2019) is the experimental gold standard:

1. Design sgRNAs tiling each candidate enhancer.
2. Transduce CRISPRi-expressing cells; FACS by gene expression (FlowFISH for endogenous; reporter for ectopic).
3. Sequence sgRNAs in low- versus high-expression bins; compute log2 enrichment per sgRNA.
4. Significance: meta-test across sgRNAs in the same enhancer.

A 2-fold expression decrease (p < 0.05) confirms the enhancer regulates the gene.

For publication-grade predictions ENCODE 4 expects:

- Test-set sensitivity and specificity against published CRISPR enhancer-screen catalogs (Fulco 2019: K562 FlowFISH; Gasperini 2019: K562; Schraivogel 2020: K562 TAP-seq). These are primarily K562 with ~thousands of pairs; the ENCODE-rE2G combined CRISPR benchmark spans ~10 cell types.
- Effect-size correlation between predicted score and observed expression effect.
- Distance-bias check (predictors over-rank close-distance pairs).

CRISPR benchmark data: the ENCODE-rE2G README points to https://github.com/EngreitzLab/CRISPR_comparison. The former `engreitzlab.org/crispri-flowfish/K562_validated_pairs.txt` URL returns an HTML page, not a table. Compute sensitivity/specificity at the calibrated threshold.

## Reconciling methods

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| ABC and ENCODE-rE2G disagree | Different feature weighting and training distributions | Both valid; report intersection as high-confidence |
| ABC strong, Cicero weak | Co-accessibility sparse for that cell type | Trust ABC if Hi-C is matched |
| HiChIP loop with no ABC prediction | Loop below ABC threshold, or peak set too narrow | Lower threshold or expand candidate enhancers |
| ENCODE-rE2G high probability, no CRISPRi support | Context-dependent biology or false positive | Prioritize for follow-up; not a publishable claim alone |
| EpiMap pair not in ABC | Reference is cell-type-aggregated | Use ABC for cell-type-specific |

## Combining predictions

`scripts/combine_predictions.py` intersects the thresholded ABC table (`EnhancerPredictionsFull_threshold*.tsv`) with the thresholded ENCODE-rE2G table on (`name`, `TargetGene`) and, with `--loops`, flags links spanned by a HiChIP loop (one anchor on the enhancer, the other on the TSS). `name` encodes coordinates, so the merge is valid only when both methods used the same candidate regions.

## Variant-to-gene use

For a non-coding GWAS lead SNP, identify the candidate enhancer (ATAC + H3K27ac peak overlap), then use ABC to find the likely target gene. Combined with deep-learning-atac (chromBPNet) variant-effect predictions, SNP at predicted enhancer + ABC target gene + chromBPNet effect is a strong functional hypothesis.

## References

- Fulco CP et al 2019 Nat Genet 51:1664 (ABC; CRISPRi-FlowFISH validation)
- Nasser J et al 2021 Nature 593:238 (ABC genome-wide application)
- Gschwind AR et al 2026 "An encyclopedia of human enhancer-gene regulatory interactions" (ENCODE-rE2G v1.0.0; preprint bioRxiv 2023.11.09.563812)
- Mumbach MR et al 2017 Nat Genet 49:1602 (HiChIP H3K27ac)
- Bhattacharyya S et al 2019 Nat Commun 10:4221 (FitHiChIP)
- Boix CA et al 2021 Nature 590:300 (EpiMap)
- Gasperini M et al 2019 Cell 176:377 (CRISPRi at scale)
- Schraivogel D et al 2020 Nat Methods 17:629 (TAP-seq enhancer screen, K562)
- Pliner HA et al 2018 Mol Cell 71:858 (Cicero)
