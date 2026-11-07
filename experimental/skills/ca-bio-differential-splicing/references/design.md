## Statistical Model Taxonomy

| Tool | Model | Test statistic | Min reps per group | Calibration regime | Fails when |
|------|-------|-----------------|---------------------|---------------------|------------|
| rMATS-turbo | Binomial counts with hierarchical PSI variance | LRT on \|ΔPSI\| > `cutoff` (default 0.0001) | n>=3 | Well-calibrated at n>=3 with adequate junction reads | Junction read imbalance; very low coverage; uncorrected for confounders |
| leafcutter | Dirichlet-multinomial GLM at cluster level | LRT on group factor | n>=3 (n=2 runs with `-i 2 -g 2`; n>=4 for calibrated calls) | 0 calls on the null at n=2-5; false calls among true ones: 4/26 at n=2, 5/30 at n=3, 1/26 at n=4 (simulation) | Undersampled clusters (DM dispersion unstable); cluster topology arbitrariness; default flags stop below 5 samples per group |
| MAJIQ deltapsi | Beta-binomial bootstrap -> posterior over PSI per LSV | P(\|ΔPSI\| > T) threshold (T=0.2) | n>=3 | Replicate-structured n=3 vs n=3 | Cohorts where between-sample variability dominates between-group |
| MAJIQ HET | Same model, heterogeneity-aware | Per-LSV permutation-based test | n>=10 | n>=10 vs n>=10 cohort designs | Tightly-controlled small replicate experiments |
| SUPPA2 (empirical) | Empirical null from between-replicate ΔPSI | ECDF on \|ΔPSI\| conditioned on TPM | n>=4 | n>=4 vs n>=4 with paired-end deep sequencing | n<=3 vs n<=3: has power but false calls inflate (see SUPPA2 section) |
| SUPPA2 (classical) | Wilcoxon rank-sum on PSI distributions | Wilcoxon p-value | n>=4 | n>=4 vs n>=4 | n<=3: smallest reachable p is 0.077 (3v3) / 0.22 (2v2), so nothing can be called |
| Shiba (2025) | Beta-binomial with explicit junction-imbalance correction | LRT | n>=2 | n=2-3 vs n=2-3 (claimed; not reproduced, see Shiba section) | Established benchmarks limited (new tool); needs XS-tagged BAMs |
| LeafcutterMD | Dirichlet-multinomial outlier mode | Per-sample p-value | n=1 vs cohort >=20 | Single-patient vs cohort | Too few controls (<20) |
| FRASER 2.0 | Beta-binomial autoencoder on Intron Jaccard Index | Per-sample p-value with delta cutoff | n=1 vs cohort >=20 | n>=20 control cohort, single-patient query | See `outlier-splicing-detection` for this regime |

MAJIQ rows are from the documentation and citations only (not run here).

The first decision is which **regime** the design falls into: between-group with replicates, heterogeneous cohort, or single-sample-vs-cohort. Within each regime, tool choice is much smaller (1-2 options).

Comprehensive 2023-2026 benchmarks: Olofsson 2023 *Biochem Biophys Res Commun*; Tran 2025 *WIREs RNA*; Kubota 2025 *NAR*. Methodology evolves — verify benchmarks and tool docs before reporting. Default 2026 recommendation: run **two complementary tools** (rMATS + leafcutter) and require concordance for high-confidence calls.

## Decision Tree by Experimental Design

| Scenario | Recommended tool | Why | Threshold |
|----------|------------------|-----|-----------|
| Standard n=3 vs n=3, GENCODE-annotated | rMATS-turbo + leafcutter (concordance) | Two algorithmic families; concordant hits = high-confidence | FDR<0.05, \|ΔPSI\|>0.10 |
| n=2 vs n=2 small pilot | Shiba or leafcutter (`-i 2 -g 2`) | Both run at n=2; Shiba's junction-imbalance correction is its claimed low-coverage advantage | FDR<0.10, \|ΔPSI\|>0.10 |
| n=10+ vs n=10+ heterogeneous (clinical, GTEx-style) | MAJIQ V3 HET | HET designed for between-sample heterogeneity | P(\|ΔPSI\|>0.2)>0.95 |
| Single rare-disease patient vs panel of n>=20 | FRASER 2.0 (see outlier-splicing-detection) | Outlier detection statistical model is fundamentally different | padj<0.05, \|delta-jaccard\|>=0.1 |
| Time-course / multi-condition design | Custom DEXSeq or limma on PSI matrix | rMATS/leafcutter primarily 2-group | FDR<0.05 on time:group interaction |
| Paired tumor-normal | rMATS with `--paired-stats` (needs PAIRADISE; see rMATS section) | Paired test reduces inter-patient variance | FDR<0.05, paired \|ΔPSI\|>0.10 |
| Cancer with spliceosomal mutation (SF3B1, U2AF1) | leafcutter or MAJIQ denovo | Cryptic events not in annotation | FDR<0.05; check 3'ss shifts in IGV |
| TDP-43 loss / ALS post-mortem | leafcutter denovo | Cryptic exons not in annotation | FDR<0.05; expect UNC13A, STMN2 |
| Non-model organism without GENCODE-grade annotation | leafcutter | Annotation-free | FDR<0.05, \|ΔPSI\|>0.10 |
| Long-read available | rMATS-long, FLAIR diffSplice | See long-read-splicing | Tool-specific |

## Design Gate (check before running anything)

- **n=1 per group: stop.** There is no between-replicate variance; rMATS still runs and returned 4 false calls (of 118 events) on a null comparison and 4 false among 22 calls on a planted one. Add replicates or use the outlier regime (`outlier-splicing-detection`).
- **Single-end reads:** rMATS needs `-t single`; the paired default gives a header-only table with exit 0.
- **Batch identical to condition** cannot be adjusted for (see Confounder Handling).
- **n<=3 per group:** not SUPPA2 (see its section); leafcutter needs `-i`/`-g` lowered.
