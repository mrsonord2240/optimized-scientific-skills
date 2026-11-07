## DGE vs DTE vs DTU: Which Question Is Being Asked?

| Question | Statistic | Tool | Example claim |
|----------|-----------|------|----------------|
| **DGE** Does the gene total change? | Sum of transcript counts | DESeq2, edgeR, limma-voom | "Gene X is upregulated 2-fold" |
| **DTE** Does this transcript change in absolute abundance? | Per-transcript count | swish (fishpond), DESeq2 on transcripts, sleuth | "Transcript X-201 is upregulated 2-fold" |
| **DTU** Do proportions of transcripts within the gene shift? | Vector of per-transcript proportions | DRIMSeq, DEXSeq, satuRn (+ stageR) | "Gene X switches from isoform 201 (50% -> 10%) to 202 (50% -> 90%)" |

DTU is statistically harder than DGE because:
1. The null is **compositional** (proportions sum to 1; one transcript up means another down).
2. **Multi-stage testing** is required: gene-level "any DTU" + transcript-level "which transcript" -> stageR formalizes this.
3. **Quantification uncertainty propagates** when transcripts are similar (Salmon EM ambiguity).

DTU and event-level differential splicing answer related but distinct questions: rMATS' `IncLevelDifference` is essentially a 1-D projection of a DTU shift onto a single event coordinate. The pragmatic 2026 default: run both an event-level tool (rMATS or leafcutter) and a DTU pipeline; reconcile.

## Tool Selection for DTU

| Tool | Model | When to use | Fails when |
|------|-------|-------------|------------|
| IsoformSwitchAnalyzeR | Wraps DEXSeq or satuRn + functional consequence annotation | Standard interpretation workflow with NMD/domain output | Manual DTU control needed; very large cohorts (>200) |
| DRIMSeq | Dirichlet-multinomial on transcript counts; gene-level DTU | Pre-filter step before DEXSeq/satuRn | Cannot annotate functional consequences alone |
| DEXSeq | Negative-binomial GLM on exon-bin or transcript counts | Classic DTU; conservative; <=5 replicates per condition | Slow at scale; uses bins not transcripts in default mode |
| satuRn | Quasi-binomial GLM with empirical-Bayes shrinkage | DTU at scale (single-cell, large bulk cohorts) | Newer; less battle-tested than DEXSeq |
| swish (fishpond) | Non-parametric SAMseq across Salmon Gibbs samples | DTE/DGE incorporating quantification uncertainty | Requires Gibbs samples; not strictly DTU |
| stageR | Two-stage testing framework | Required for proper OFDR control on top of DRIMSeq/DEXSeq/satuRn | Standalone — wraps another tool's output |
| sleuth | Bootstrap-based DTE on kallisto | When committed to kallisto pipeline | Less active development; superseded by fishpond+swish |

`isoformSwitchAnalysisPart1()` applies the <=5 / >5 rule automatically. This Skill calls the direct wrappers deliberately: use `isoformSwitchTestDEXSeq()` when every condition has <=5 replicates, and `isoformSwitchTestSatuRn()` when any has >5. The direct 2.6.0 wrappers only warn on a mismatch (DEXSeq: "You seem to have many replicates... use isoformSwitchTestSatuRn()"; satuRn: "You seem to have few replicates... use isoformSwitchTestDEXSeq()"). At exactly 5 replicates results may differ; document the choice.

## Decision Tree by Research Question

| Question | Recommended approach |
|----------|----------------------|
| Functional consequences of switches (domains, NMD, signal peptide) | IsoformSwitchAnalyzeR with the annotators available to you (CPC2 + Pfam locally; SignalP/IUPred2A/DeepTMHMM if licensed) |
| Pure statistical DTU (gene-level + transcript-level OFDR) | DRIMSeq (filter) -> DEXSeq -> stageR; or -> satuRn -> stageR for n>5 |
| DTU with proper quantification uncertainty | Salmon `--numGibbsSamples 20` -> tximeta -> swish for DTE; concurrent DTU |
| Single-cell DTU | satuRn (DEXSeq doesn't scale to scRNA-seq) |
| Long-read DTU (PacBio Iso-Seq, ONT) | IsoformSwitchAnalyzeR on the long-read count matrix via `importRdata` (no Salmon EM uncertainty) |
| Time-course DTU | DEXSeq with time as factor + interaction; or limma::lmFit on logit-prop matrix |
| Cancer / disease — switch hits -> mechanism | Standard pipeline + cross-reference with eCLIP, ClinVar, COSMIC |
| Therapeutic ASO target identification | Standard pipeline + sashimi visualization + SpliceAI design |
