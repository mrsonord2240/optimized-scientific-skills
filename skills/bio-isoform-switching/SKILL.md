---
name: bio-isoform-switching
description: Analyzes differential transcript usage (DTU) and isoform switches with functional consequence prediction (NMD via 50nt rule, ORF disruption, protein domain loss/gain, coding-potential shifts; signal peptide and IDR changes when licensed annotators are available). Tools include IsoformSwitchAnalyzeR (DEXSeq up to 5 replicates per condition, satuRn above), the manual DRIMSeq -> DEXSeq/satuRn -> stageR DTU pipeline, and fishpond/swish for inferential-uncertainty-aware DTE. Distinguishes DTU from DGE and DTE; integrates external annotators (CPC2 and Pfam run locally; SignalP, IUPred2A and DeepTMHMM are licence-gated). Use when investigating how splicing differences alter protein function or trigger NMD-mediated degradation.
tool_type: r
primary_tool: IsoformSwitchAnalyzeR
license: MIT
author: GPTomics
---

## Version Compatibility

Checked 2026-09 on R 4.4.3 / Bioconductor 3.20: IsoformSwitchAnalyzeR 2.6.0, DRIMSeq 1.34.0, DEXSeq 1.52.0, satuRn 1.14.0, stageR 1.28.0, fishpond 2.12.0, tximport 1.34.0, tximeta 1.24.0, on Salmon 2.7.0 output. CPC2 (standalone), HMMER 3.4 with Pfam-A (2026-09 release). The IsoformSwitchAnalyzeR v2 paper (Han et al. 2026, *NAR Genom Bioinform*; preprint 2025) says v2 uses DEXSeq for smaller studies and satuRn for larger ones and for single-cell data, and imports quantifications from long-read and single-cell pipelines. **2.6.0 does not choose the test for you** (it only warns; call the DEXSeq or satuRn test explicitly, see below), and no newer release was available to check whether a later one does.

```r
BiocManager::install(c('IsoformSwitchAnalyzeR', 'DRIMSeq', 'DEXSeq', 'satuRn', 'stageR', 'fishpond', 'tximeta', 'tximport'))
```
```bash
conda install -c bioconda salmon hmmer        # Salmon with --numGibbsSamples 20 for swish; hmmscan for Pfam
# CPC2 (standalone, Python 2 code), SignalP / IUPred2A / DeepTMHMM (licence-gated): see "Functional Consequence Annotation"
```

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Isoform Switching and Differential Transcript Usage

Identify shifts in *which* transcript a gene predominantly uses between conditions, and predict functional consequences. Statistically distinct from DGE and DTE; biologically distinct because the same gene-level expression can hide a complete isoform switch with major protein-level consequences.

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

**Choose the test yourself: `isoformSwitchTestDEXSeq()` when every condition has <=5 replicates, `isoformSwitchTestSatuRn()` when any has >5.** IsoformSwitchAnalyzeR 2.6.0 enforces this only by warning (DEXSeq wrapper: "You seem to have many replicates... use isoformSwitchTestSatuRn()"; satuRn wrapper: "You seem to have few replicates... use isoformSwitchTestDEXSeq()"). At exactly 5 replicates results may differ; document the choice.

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

## IsoformSwitchAnalyzeR Workflow

**Goal:** Identify isoform switches with functional consequences in one integrated workflow.

**Approach:** Import Salmon (raw counts), join the design to the quantification by sample name, pre-filter, run the statistical test, then annotate consequences (next section).

```r
library(IsoformSwitchAnalyzeR)

salmonQuant <- importIsoformExpression(
    parentDir = 'salmon_quant/',              # one sub-directory per sample, each holding quant.sf
    calculateCountsFromAbundance = FALSE,     # raw NumReads as counts; see "Count route"
    addIsofomIdAsColumn = TRUE                # (sic) the package's own spelling
)

# Join the design to the quantification BY SAMPLE NAME. importIsoformExpression sorts samples
# alphabetically, so a hand-typed condition vector silently mislabels samples (on SRR-style IDs it
# recovered 0/20 planted switches with no warning; joined by name, 20/20).
meta <- read.delim('sample_metadata.tsv')     # columns: sample_id (= quant sub-directory name), condition, [batch]
ids <- setdiff(colnames(salmonQuant$counts), 'isoform_id')
if (!setequal(ids, meta$sample_id))
    stop('sample IDs differ between quantification and metadata.\n  only in quant: ', paste(setdiff(ids, meta$sample_id), collapse = ', '),
         '\n  only in metadata: ', paste(setdiff(meta$sample_id, ids), collapse = ', '))
stopifnot('duplicated sample_id in metadata' = !anyDuplicated(meta$sample_id))
design <- data.frame(sampleID = ids, condition = meta$condition[match(ids, meta$sample_id)])
# design$batch <- meta$batch[match(ids, meta$sample_id)]   # extra columns = covariates, see below
print(table(design$condition))

aSwitchList <- importRdata(
    isoformCountMatrix = salmonQuant$counts,
    isoformRepExpression = salmonQuant$abundance,
    designMatrix = design,
    isoformExonAnnoation = 'annotation.gtf',   # (sic); a GTF with CDS lines also supplies annotated ORFs
    isoformNtFasta = 'transcripts.fa',
    addAnnotatedORFs = TRUE,
    showProgress = FALSE
)

aSwitchList <- preFilter(
    aSwitchList,
    geneExpressionCutoff = 1,
    isoformExpressionCutoff = 0,
    IFcutoff = 0.01,
    removeSingleIsoformGenes = TRUE,
    keepIsoformInAllConditions = TRUE
)

# DEXSeq for <=5 replicates per condition, satuRn above (see Tool Selection).
# reduceToSwitchingGenes = FALSE: with TRUE (the function default) a run with no switches stops with an error
aSwitchList <- if (max(table(design$condition)) > 5) {
    isoformSwitchTestSatuRn(aSwitchList, reduceToSwitchingGenes = FALSE, alpha = 0.05, dIFcutoff = 0.1, diagplots = FALSE)
} else {
    isoformSwitchTestDEXSeq(aSwitchList, reduceToSwitchingGenes = FALSE, alpha = 0.05, dIFcutoff = 0.1)
}
f <- aSwitchList$isoformFeatures
sum(f$isoform_switch_q_value < 0.05 & abs(f$dIF) > 0.1, na.rm = TRUE)   # 0 -> see "Count route"
```

`examples/isoform_switch_analysis.R` runs this whole workflow, including consequences and the plot, on a built-in synthetic dataset (`Rscript isoform_switch_analysis.R`) or on your own Salmon directory, GTF, FASTA and sample-metadata table.

`preFilter` parameters:
- `geneExpressionCutoff = 1` — minimum TPM for gene to be tested (raise for stricter)
- `isoformExpressionCutoff = 0` — minimum TPM per isoform (set to 1 for stricter)
- `IFcutoff = 0.01` — minimum isoform fraction; below = noise
- `removeSingleIsoformGenes = TRUE` — drop genes with only one detectable isoform (cannot have DTU)
- `keepIsoformInAllConditions = TRUE` — require expression across all conditions

**Replicates and covariates.** IsoformSwitchAnalyzeR stops on 1 replicate per condition ("A statistical test cannot be performed without replicates"). 2 per condition runs but is not trustworthy on the raw-count route (see "Count route": label-shuffled 2 v 2 splits of real data are called almost as often as the true split); plan >=3, ideally >=5. Any extra `designMatrix` column is treated as a covariate: on a batch-confounded synthetic set a batch column removed the batch artefacts (19/30 artefact genes called without it, 0 with it) at no cost to the 20/20 planted switches; a covariate identical to condition stops with "not full rank", and a constant one with "Contain constant information". `importRdata(detectUnwantedEffects = TRUE)` (the default) also runs `sva` and adds any surrogate variables it finds to the model (`sv1` appeared for one of the label-shuffled chrX 2 v 2 splits; check `colnames(aSwitchList$designMatrix)`). Per `?importRdata` and the package NEWS, abundance, IF and dIF are then corrected for these covariates (counts stay unmodified but the test includes them), so with a covariate the reported dIF is no longer the plain mean of per-sample isoform fractions (up to 0.08 from a hand computation on the synthetic 3 v 3 set with a batch column; 7e-5 without).

For long-read data pass the long-read transcript count matrix to `importRdata` directly (a count-only import recovered 20/20 planted switches on the synthetic set) — no Salmon EM uncertainty.

### Count route (`calculateCountsFromAbundance`)

`importIsoformExpression()` derives counts from abundance by default (`calculateCountsFromAbundance = TRUE` = tximport `scaledTPM`); the package vignette recommends this because it carries Salmon's bias correction into the counts. The tests use those counts, while dIF comes from the abundance in either route. **On real chrX data this default called nothing in every label split (below), so it cannot be validated there. This Skill uses raw NumReads (`FALSE`), which does find switches but is not validated either: label-shuffled splits of the same data are called almost as often as the true one. Read "Null check" before trusting a call.**

Real chrX RNA-seq, 2 GBR v 2 YRI (nf-core rnasplice, Salmon 2.7.0, GRCh37), `isoformSwitchTestDEXSeq`, switching isoform = q < 0.05 and |dIF| > 0.1:

| Count route | Switching isoforms / genes | RPL10 ENST00000406022 (dIF -0.33) q |
|-------------|----------------------------|-------------------------------------|
| `calculateCountsFromAbundance = FALSE` (raw NumReads) | 26 / 20 | 4.6e-18 |
| tximport `lengthScaledTPM` | 28 / 22 | 1.9e-18 |
| default (`scaledTPM`) | 0 / 0 | 1 |
| tximport `dtuScaledTPM` | 0 / 0 | 1 |

- **Cause (consistent with, not proven by, the runs below):** scaledTPM counts are proportional to TPM, i.e. each isoform's counts are divided by its length. The factor median(length)/length ranges 0.019-576 across the 5,895 isoforms with `EffectiveLength` (0.023-18 with `Length`), so short isoforms are inflated and stop behaving like counts, and DEXSeq's isoform-versus-rest test degrades. Multiplying the raw counts by either factor (`Length` or `EffectiveLength`) reproduced the result (0 switches, RPL10 q = 1); multiplying by the same factors shuffled across isoforms did not (34 switches / 27 genes with `Length`, 31 / 23 with `EffectiveLength`). satuRn was not compared across routes on real data; the manual-pipeline route numbers are under "Manual DTU Pipeline".
- The synthetic planted set has near-uniform isoform lengths and recovered 20/20 with either route, so a clean simulation will not reveal the problem.

### Null check: shuffle the condition labels

The four chrX samples are 2 GBR + 2 YRI. Splits that put one GBR and one YRI in each group carry no population signal, so any calls there are individual-level differences that the n = 2 dispersion estimate does not absorb. Genes with a switching isoform (q < 0.05, |dIF| > 0.1), raw counts unless stated:

| Split | ISAR DEXSeq | ISAR default route | manual DEXSeq gene q < 0.05 | manual DRIMSeq | manual both | ISAR calls also DEXSeq-confirmed |
|-------|-------------|--------------------|-----------------------------|----------------|-------------|----------------------------------|
| true GBR v YRI | 20 | 0 | 13 | 11 | 7 | 12 |
| mixed 1 | 11 | 0 | 0 | 4 | 0 | 0 |
| mixed 2 | 15 | 0 | 4 | 3 | 1 | 4 |

- **The raw route is not validated at n = 2.** The mixed splits are called at 55% and 75% of the true split's count, and only 1 of the 20 true-split genes recurs in either. No cutoff separates them cleanly (genes, true / mixed 1 / mixed 2: q < 0.01: 12 / 3 / 7; q < 1e-4: 8 / 2 / 5; q < 1e-6: 5 / 2 / 1; |dIF| > 0.3: 17 / 7 / 10; a higher gene-TPM or isoform-TPM filter changed nothing). Report a 2 v 2 list as exploratory. The default route's zeros show only that it has no power on this data, not that it is the safe choice.
- **What reduced it here:** keep only calls that the manual gene-level DEXSeq test above also makes (12 of 20 true-split genes, 0 and 4 of the mixed-split calls; DRIMSeq as well: 6, 0, 1). It cuts the shuffled calls but does not remove them (4 of 15 remain), and the counts are from one 4-sample data set.
- **Run the permutation check on every real data set** (this block reuses `salmonQuant` and `design` from the workflow; each call repeats the import, about 20-40 s on chrX). Trust an observed list only when it is far above every permuted count. Group sizes of 2 give only two alternative splits (no p-value); 3 v 3 gives nine, 6 v 6 hundreds. Splits that share half of the true split keep part of the signal, so they overstate the null (conservative).

```r
# Same import, filter and test as the workflow, on any design; returns the genes called
call_genes <- function(design) {
    sl <- importRdata(isoformCountMatrix = salmonQuant$counts, isoformRepExpression = salmonQuant$abundance,
                      designMatrix = design, isoformExonAnnoation = 'annotation.gtf', isoformNtFasta = 'transcripts.fa',
                      addAnnotatedORFs = FALSE, showProgress = FALSE, quiet = TRUE)
    sl <- preFilter(sl, geneExpressionCutoff = 1, isoformExpressionCutoff = 0, IFcutoff = 0.01,
                    removeSingleIsoformGenes = TRUE, keepIsoformInAllConditions = TRUE, quiet = TRUE)
    sl <- if (max(table(design$condition)) > 5) {
        isoformSwitchTestSatuRn(sl, reduceToSwitchingGenes = FALSE, alpha = 0.05, dIFcutoff = 0.1, diagplots = FALSE, quiet = TRUE)
    } else {
        isoformSwitchTestDEXSeq(sl, reduceToSwitchingGenes = FALSE, alpha = 0.05, dIFcutoff = 0.1, quiet = TRUE)
    }
    f <- sl$isoformFeatures
    unique(f$gene_id[!is.na(f$isoform_switch_q_value) & f$isoform_switch_q_value < 0.05 & abs(f$dIF) > 0.1])
}

# Shuffled labels that keep the group sizes and share only the chance overlap with the true split
lab <- design$condition; n <- length(lab); in1 <- lab == lab[1]; k <- sum(in1)
canon <- function(v) paste(v == v[1], collapse = '')          # a split and its mirror image are the same split
seen <- canon(in1); perms <- list(); set.seed(1)
for (i in 1:5000) {
    if (length(perms) >= 10) break
    p <- sample(lab); key <- canon(p == lab[1])
    if (abs(sum(p == lab[1] & in1) - k * k / n) <= 0.5 && !key %in% seen) { seen <- c(seen, key); perms[[length(perms) + 1]] <- p }
}
observed <- length(call_genes(design))
permuted <- sapply(perms, function(p) { d <- design; d$condition <- p; tryCatch(length(call_genes(d)), error = function(e) NA) })
cat('observed', observed, 'genes; label-permuted:', permuted, '\n')
```

Checked with this block: real chrX 2 v 2 observed 20, permuted 15 and 11. Synthetic 6 v 6 (satuRn branch, batch covariate kept): observed 24 (20/20 planted, 0 null-type genes), permuted 0 0 0 1 0 0 0 0 0 0 over 10 splits. Synthetic 3 v 3 (DEXSeq): observed 25 (20/20 planted, 2 null-type genes), permuted 0 2 0 0 0 0 0 0 0 over the 9 other splits.

**Recommended route:** raw NumReads (`FALSE`); at least 3 replicates per condition, ideally 5 or more; run the permutation check; confirm calls with the manual gene-level test; if you also try `TRUE`, compare the two switch lists and inspect the top switches (`switchPlot`). A large gap between routes means the count route, not biology, drove the result.

## Functional Consequence Annotation

**Goal:** Predict how each switch alters protein structure, function, and stability.

**Approach:** Get ORFs, extract sequences, run the external annotators outside R, import their results, then request only the consequence types whose annotation was imported.

Order (each step needs the previous one; verified on a synthetic GTF without CDS and on the real chrX GTF with CDS):

1. **ORFs first.** `importRdata(addAnnotatedORFs = TRUE)` loads CDS from a GTF that has CDS lines (`orf_origin` = Annotation). Only if the GTF has no CDS (or isoforms are novel) run `analyzeORF(aSwitchList, orfMethod = 'longest', genomeObject = NULL)`. `extractSequence()` and `analyzePFAM()` stop without ORFs ("Please run the 'addORFfromGTF()' ... function(s) to detect ORFs").
2. **Do not run `analyzeORF()` over annotated ORFs** (for novel isoforms the package vignette adds `addORFfromGTF()` + `analyzeNovelIsoformORF()`; not run here). It overwrites all of them (real chrX: `orf_origin` Annotation 1738 -> Predicted 1738) and worsens the NMD call: PTC flag versus Ensembl `nonsense_mediated_decay` biotype was sensitivity 0.96 / specificity 1.00 with annotated ORFs and 0.64 / 0.95 after `analyzeORF('longest')`.
3. `analyzeAlternativeSplicing(aSwitchList, onlySwitchingGenes = TRUE)` — needed for `intron_retention`.
4. `extractSequence(aSwitchList, onlySwitchingGenes = TRUE, pathToOutput = 'sequences/', writeToFile = TRUE)` (create the directory first) writes `isoformSwitchAnalyzeR_isoform_nt.fasta` and `..._AA.fasta`.
5. Run the annotators (table) and import each result.
6. `analyzeSwitchConsequences()` with only the types you imported.

```r
# ORFs: annotated CDS come with importRdata; predict only when the GTF had none
if (!any(aSwitchList$orfAnalysis$orf_origin == 'Annotation', na.rm = TRUE)) {
    aSwitchList <- analyzeORF(aSwitchList, orfMethod = 'longest', genomeObject = NULL)
}
aSwitchList <- analyzeAlternativeSplicing(aSwitchList, onlySwitchingGenes = TRUE)
dir.create('sequences', showWarnings = FALSE)
aSwitchList <- extractSequence(aSwitchList, onlySwitchingGenes = TRUE, pathToOutput = 'sequences/', writeToFile = TRUE)

# after running the external tools on sequences/:
aSwitchList <- analyzeCPC2(aSwitchList, pathToCPC2resultFile = 'cpc2_result.txt', removeNoncodinORFs = FALSE)
aSwitchList <- analyzePFAM(aSwitchList, pathToPFAMresultFile = 'pfam_scanfmt.txt')
# licence-gated, not run here: analyzeSignalP(pathToSignalPresultFile=), analyzeIUPred2A(pathToIUPred2AresultFile=)

aSwitchList <- analyzeSwitchConsequences(
    aSwitchList,
    consequencesToAnalyze = c('intron_retention', 'ORF_seq_similarity', 'NMD_status',
                              'coding_potential', 'domains_identified'),   # + one type per extra annotator imported
    dIFcutoff = 0.1
)
```

| External tool | Purpose | Consequence types it enables | Status |
|----------------|---------|------------------------------|--------|
| CPC2 (standalone, Python 2 code) | Coding vs non-coding | `coding_potential` | Ran on `isoformSwitchAnalyzeR_isoform_nt.fasta`: `python CPC2.py -i <nt.fasta> -o cpc2_result` (needs a Python 2.7 env and its bundled libsvm built) |
| Pfam: `hmmscan --cut_ga --domtblout` against Pfam-A, then `examples/hmmscan_to_pfamscan.py` | Protein domains | `domains_identified` | Ran (HMMER 3.4 output; below) |
| SignalP 6.0 | Signal peptides | `signal_peptide_identified` | Licence-gated, **not run**; import per `?analyzeSignalP` |
| IUPred2A or NetSurfP-2 | Intrinsically disordered regions | `IDR_identified`; `IDR_type` needs IUPred2A | Licence/registration-gated, **not run**; read from source: `IDR_identified` needs `analyzeIUPred2A()` or `analyzeNetSurfP2()` |
| DeepTMHMM | Transmembrane topology | `isoform_topology` (via `analyzeDeepTMHMM()`), not the IDR types | Cloud service, **not run**; read from source only |

Behaviours checked by running:
- **Pfam import.** `analyzePFAM()` rejects a raw `hmmscan --domtblout` file ("more columns than column names") and accepts only a pfam_scan-style table with a clan accession (CL...) in column 15. Convert with the shipped script, which adds clans from Pfam-A.clans.tsv:

  ```bash
  hmmscan --cut_ga --domtblout pfam_domtbl.txt Pfam-A.hmm sequences/isoformSwitchAnalyzeR_isoform_AA.fasta > /dev/null
  curl -O https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.clans.tsv.gz
  python examples/hmmscan_to_pfamscan.py pfam_domtbl.txt pfam_scanfmt.txt Pfam-A.clans.tsv.gz
  ```
  On the synthetic set (10 genes carrying a planted ubiquitin domain in exon E3): PF00240 found in every isoform containing it and never in the exon-skipped one; 60 domain rows imported; "Domain loss" called 10/10 in the right direction.
- **Missing annotators are an error, not silent.** Requesting a type without its annotation stops: "To test differences in signal peptides, the result of the SignalP analysis must be advailable (sic)"; likewise for CPC2/CPAT (`coding_potential`), Pfam (`domains_identified`) and NetSurfP2/IUPred2A (`IDR_*`). Build `consequencesToAnalyze` only from what you imported. A missing file gives "The file(s) ... does not exist".
- **`removeNoncodinORFs` has no default.** `TRUE` drops the ORFs of CPC2-"noncoding" isoforms, which removes PTC-bearing isoforms from the NMD call: on the planted set 8/10 poison isoforms were called NMD-sensitive with `TRUE` (2 were CPC2 noncoding) and 10/10 with `FALSE` or no CPC2. Use `FALSE` when NMD is the question.

The external tools run *outside* R; IsoformSwitchAnalyzeR writes the FASTA files and re-imports the parsed results.

## NMD Prediction (The 50-nt Rule)

A transcript is predicted NMD-sensitive if its premature termination codon (PTC) lies **>50-55 nt upstream of the last exon-exon junction** (Maquat 2004 *Nat Rev Mol Cell Biol*; Lykke-Andersen & Jensen 2015 *Nat Rev Mol Cell Biol*).

**Mechanism:** Spliceosome deposits the Exon Junction Complex (EJC) ~20-24 nt upstream of every exon-exon junction. During the pioneer round of translation, ribosome reading through removes EJCs upstream of the stop codon. If a stop codon precedes the last EJC by >50 nt, the EJC remains, recruits UPF1 -> SMG1 phosphorylation -> SMG6/SMG7 -> mRNA decay.

**Caveats and exceptions:**
- **Last-exon PTCs escape NMD** — can be dominant-negative or gain-of-function (e.g. MYH7 truncating variants).
- **3'UTR length matters**: very long 3' UTRs (>1 kb past stop) trigger NMD via UPF1 binding even without EJCs (faux-3'UTR rule).
- **Tissue-specific NMD**: SMG6 vs SMG5/7 ratios vary; UPF1 stress conditions modulate.
- **PTC distance must be measured on the spliced transcript**, not the genomic distance.
- **Some predicted-NMD transcripts escape decay in RNA-seq** (Lindeboom 2016 *Nat Genet*; the escape fraction is tissue- and rule-dependent and was not verified here). Treat NMD prediction as probabilistic, not certain.

`analyzeSwitchConsequences` with `'NMD_status'` evaluates this from the ORF (annotated if the GTF has CDS, else predicted; see the ORF order above) and the transcript model. Checked: in the shipped demo the PTC flag marks exactly the 8 planted poison isoforms, and on real chrX the annotated-ORF PTC flag agrees with Ensembl's NMD biotype (sensitivity 0.96, above).

## AS-NMD as a Regulatory Layer

A large class of conserved alternative splicing events is **deliberately PTC-introducing** to titrate functional protein levels:

- **All major SR proteins** (SRSF1-12) autoregulate via poison exons (Lareau 2007 *Nature*; Ni 2007 *Genes Dev*)
- **All major hnRNPs** likewise
- **Ribosomal protein genes** use AS-NMD autoregulation (e.g. rpL3, rpL12; Cuccurese 2005 *NAR*)
- **SCN1A** poison exon -> STK-001 (Stoke Therapeutics) ASO for Dravet syndrome, in clinical development (Han 2020 *Sci Transl Med*; trial phase not verified here)

**Functional implication:** an *increase* in PSI of a poison exon *decreases* functional protein. Sign-of-effect in DTU output is opposite from intuition for these genes. Always check whether the alternative form is PTC-bearing before interpreting direction.

**Disease examples:**
- TDP-43 cryptic exon in UNC13A introduces a PTC -> NMD of a disease-relevant transcript (Brown 2022 *Nature*). The STMN2 cryptic exon works differently: premature polyadenylation gives a truncated transcript, not PTC-NMD (Melamed 2019 *Nat Neurosci*; Klim 2019 *Nat Neurosci*). These appear in post-mortem ALS brain, not blood: check the tissue.
- Last-exon truncating variants in TTN (and MYH7): escape NMD -> stable poison/dominant-negative protein

## Manual DTU Pipeline (DRIMSeq + DEXSeq + stageR)

The canonical reference is the *F1000Research* "Swimming downstream" workflow (Love, Soneson, Patro 2018; Bioconductor `rnaseqDTU`). Below: tximport with raw counts (`countsFromAbundance = 'no'`, the route that worked on real data; see Count route).

```r
library(tximport); library(DRIMSeq); library(DEXSeq); library(stageR)

# meta: data.frame(sample_id, condition); files: named vector of quant.sf paths (names = sample_id)
# tx2gene: data.frame(tx = transcript ID, gene = gene ID)
txi <- tximport(files, type = 'salmon', txOut = TRUE, countsFromAbundance = 'no')
cts <- txi$counts
txdf <- data.frame(gene_id = tx2gene$gene[match(rownames(cts), tx2gene$tx)], feature_id = rownames(cts), cts, check.names = FALSE)
txdf <- txdf[!is.na(txdf$gene_id), ]

samples <- data.frame(sample_id = colnames(cts), condition = factor(meta$condition[match(colnames(cts), meta$sample_id)]))
stopifnot(!anyNA(samples$condition))

n <- nrow(samples); n_small <- min(table(samples$condition))
d <- dmDSdata(counts = txdf, samples = samples)
d <- dmFilter(d, min_samps_gene_expr = n, min_gene_expr = 10,
              min_samps_feature_expr = n_small, min_feature_expr = 10,
              min_samps_feature_prop = n_small, min_feature_prop = 0.1)

# Qualify with DRIMSeq:: -- once DEXSeq is attached, Biobase::samples masks DRIMSeq::samples() ("unable to find an
# inherited method for function 'samples' for signature 'object = "dmDSdata"'"); counts() is also exported by
# DEXSeq, DESeq2 and BiocGenerics, so qualify it too.
cnt <- DRIMSeq::counts(d)
dxd <- DEXSeqDataSet(
    countData = round(as.matrix(cnt[, samples$sample_id])),
    sampleData = DRIMSeq::samples(d),
    design = ~ sample + exon + condition:exon,
    featureID = cnt$feature_id,
    groupID = cnt$gene_id
)
dxd <- estimateSizeFactors(dxd)
dxd <- estimateDispersions(dxd, quiet = TRUE)
dxd <- testForDEU(dxd, reducedModel = ~ sample + exon)
dxr <- DEXSeqResults(dxd, independentFiltering = FALSE)
qval <- perGeneQValue(dxr)

pConfirmation <- matrix(dxr$pvalue, ncol = 1, dimnames = list(dxr$featureID, NULL))
pConfirmation[is.na(pConfirmation)] <- 1
tx2gene_d <- data.frame(transcript = dxr$featureID, gene = dxr$groupID)

stageRObj <- stageRTx(pScreen = qval, pConfirmation = pConfirmation, pScreenAdjusted = TRUE, tx2gene = tx2gene_d)
stageRObj <- stageWiseAdjustment(stageRObj, method = 'dtu', alpha = 0.05)

results <- getAdjustedPValues(stageRObj, order = FALSE, onlySignificantGenes = FALSE)   # geneID, txID, gene, transcript
```

Confirm IsoformSwitchAnalyzeR calls (`aSwitchList` from the workflow) with the gene-level test above; see "Count route" for why:

```r
f <- aSwitchList$isoformFeatures
hit <- unique(f$isoform_id[!is.na(f$isoform_switch_q_value) & f$isoform_switch_q_value < 0.05 & abs(f$dIF) > 0.1])
called <- unique(tx2gene$gene[match(hit, tx2gene$tx)])            # gene IDs of the tx2gene table, not gene names
confirmed <- intersect(called, names(qval)[!is.na(qval) & qval < 0.05])
```

Checked (synthetic 6 v 6, planted truth): `perGeneQValue` < 0.05 in 27 genes = 20/20 planted + 1 null-gene false positive + 5 planted 0.15-shift and 1 planted 0.06-shift genes; DRIMSeq (`dmPrecision`/`dmFit`/`dmTest`) called 27 genes, Jaccard 0.93 with DEXSeq; stageR confirmed the truly switching transcript (poison B or skip C) in 20/20. Control-versus-control split: 2 gene calls of 291. Real chrX 2 v 2 (same filter; raw / `lengthScaledTPM` / `scaledTPM` counts): DEXSeq 13 / 13 / 2 genes; DRIMSeq 11 / 11 / 7; DEXSeq and DRIMSeq overlap 7 genes with raw counts. Label-shuffled splits: see "Null check".

**Importing with tximeta instead.** `tximeta(coldata)` needs a linked transcriptome that matches the Salmon index; on Salmon 2.7.0 with tximeta 1.24.0 it did not match ("couldn't find matching transcriptome") and returned an empty `rowData`, so `rowData(se)$gene_id` / `tx_id` could **not** be run here. Simulated only: `gene_id` is a list-like column, so `unlist()` it before building the `data.frame` above. `tximeta(coldata, skipMeta = TRUE)` works (used for swish below) but carries no gene mapping.

**stageR semantics:**
- **Stage 1 (screening)**: gene-level p-value (`perGeneQValue` from DEXSeq, or DRIMSeq's gene-level adjusted p) is filtered at the desired Overall FDR.
- **Stage 2 (confirmation)**: only within significant genes, individual transcripts are tested at a within-gene FWER computed to maintain global OFDR.
- **Net effect**: gene-level FDR is properly controlled, AND the transcript that drove the call is known.
- Without stageR: naive transcript-level BH overcounts because the gene-level multiple-testing burden is ignored.

## fishpond/swish for Inferential-Uncertainty-Aware Testing

**Goal:** Test DTE while propagating quantification uncertainty from Salmon's Gibbs samples.

**Approach:** Run Salmon with `--numGibbsSamples 20`, import with tximeta, then use swish to average a non-parametric SAMseq-style test across inferential replicates.

```r
library(fishpond); library(tximeta)

# meta: data.frame(sample_id, condition); salmon_dir/<sample_id>/quant.sf from a run with --numGibbsSamples 20
coldata <- data.frame(names = meta$sample_id,
                      files = file.path(salmon_dir, meta$sample_id, 'quant.sf'),
                      condition = factor(meta$condition))   # swish needs a factor
se <- tximeta(coldata, skipMeta = TRUE)    # drop skipMeta once a linkedTxome matches your index
y <- scaleInfReps(se)
y <- labelKeep(y)
y <- y[mcols(y)$keep, ]

set.seed(1)
y <- swish(y, x = 'condition')
y <- computeInfRV(y)                        # adds mcols(y)$meanInfRV

dte_results <- as.data.frame(mcols(y))     # log2FC, pvalue, qvalue, meanInfRV
sig <- subset(dte_results, qvalue < 0.05)
```

Checked on real chrX Gibbs samples (fishpond 2.12.0, 2 v 2): 391 transcripts tested, 6 at q < 0.05 (top ENST00000380861, log2FC -1.45 with NumReads 75/63 vs 9/34); with a character `condition` swish stops with "is.factor(condition) is not TRUE". On a synthetic 6 v 6 set with planted DTE, swish called 33/35 planted up-transcripts (all log2FC > 0) with 29/848 false positives.

**`infRV`** (inferential relative variance) is a per-feature uncertainty diagnostic (`meanInfRV`); high-infRV transcripts are unreliable and can be filtered before testing (e.g. `y[mcols(y)$meanInfRV < 1, ]`, choose the cut from `summary(mcols(y)$meanInfRV)`; median 0.69 on the chrX set). Critical for genes with many similar isoforms (TTN, MAPT, NEFM) where Salmon's EM is uncertain.

## Per-Tool Failure Modes

### DEXSeq: Slowness at Scale

**Trigger:** Bulk cohort with >50 samples or single-cell DTU.

**Mechanism:** DEXSeq fits a NB GLM per exon-bin per gene; computational cost scales linearly with samples × bins.

**Symptom:** `estimateDispersions` takes hours; `testForDEU` exhausts memory.

**Fix:** Switch to satuRn (designed for scale, including scRNA-seq); run with parallelization (`BPPARAM = MulticoreParam(8)`).

### DRIMSeq: Filtering Sensitivity

**Trigger:** `dmFilter` thresholds that do not fit the cohort.

**Mechanism:** `dmFilter`'s own defaults are all 0 (no filtering; DRIMSeq 1.34.0). The values in the manual pipeline are suggestions scaled to group size: features need >= 10 counts in `n_small` samples and proportion >= 0.1 in `n_small` samples. `min_samps_*` larger than the number of samples stops with "min_samps_gene_expr >= 0 && min_samps_gene_expr <= ncol(x@counts) is not TRUE", and thresholds nothing passes with "!No genes left after filtering!".

**Symptom:** Most candidate genes filtered out; few testable genes.

**Fix:** Tune to dataset: lower thresholds for low-coverage data, raise for high-coverage. Document choice.

### satuRn: Empirical-Bayes Shrinkage Limits

**Trigger:** Very small cohort (n=2 vs n=2) or very heterogeneous.

**Mechanism:** Empirical-Bayes shrinkage assumes shared dispersion across genes; collapses with too-few or too-heterogeneous samples.

**Symptom:** Inflated p-values; few discoveries despite real effects.

**Fix:** Aggregate replicates (pseudobulk), or switch to DEXSeq for small cohorts; use larger cohorts when possible.

### swish: Salmon Gibbs Requirements

**Trigger:** Running swish on Salmon output without Gibbs samples.

**Mechanism:** swish averages over inferential replicates from Salmon's Gibbs sampler; requires `--numGibbsSamples 20` (or bootstrap with `--numBootstraps`) at Salmon time.

**Symptom:** `scaleInfReps` stops with "there are no inferential replicates in the assays of 'y'".

**Fix:** Re-run Salmon with `--numGibbsSamples 20`; this triples Salmon runtime but enables uncertainty-aware testing.

## Reconciliation: When DTU and Event-Level Disagree

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| Significant DTU, no rMATS hit | DTU shift across many transcripts; no single canonical event captures it | Examine isoform structure in switchPlot; report at gene level |
| rMATS sig, no significant DTU | Single event in single isoform; not a gene-level DTU | Report as event-level result; DTU not the right framing |
| Both sig, same gene, different "main" isoforms | Annotation differs (rMATS uses GENCODE basic; ISA uses comprehensive) | Standardize annotation; re-run |
| DTU shows poison-exon switch, gene-level DGE shows decrease | NMD-coupled regulation: AS-NMD reducing protein on top of transcription | Mechanism: AS-NMD; report direction carefully |

For high-confidence reporting: concordant DTU + event-level + sashimi visualization.

## Single-Cell DTU

For scRNA-seq, **satuRn** scales where DEXSeq does not (satuRn has explicit single-cell calibration, Gilis 2022 *F1000Research*). IsoformSwitchAnalyzeR 2.6.0's `importRdata` has no single-cell argument; a single-cell count matrix route was not tested here.

**Strong recommendation:** pseudobulk by cell type first; per-cell DTU is rarely powered with droplet 3' chemistry. See `single-cell-splicing` for chemistry-specific limitations.

## Visualization

```r
extractTopSwitches(
    aSwitchList,
    filterForConsequences = TRUE,
    n = 25,
    sortByQvals = TRUE
)

switchPlot(
    aSwitchList,
    gene = 'TARGET_GENE',
    condition1 = 'control',
    condition2 = 'treatment',
    localTheme = theme_bw(base_size = 12)
)

extractSwitchSummary(aSwitchList, filterForConsequences = TRUE)
extractConsequenceSummary(aSwitchList, consequencesToAnalyze = 'all', plotGenes = FALSE)
extractConsequenceEnrichment(aSwitchList, consequencesToAnalyze = 'all')
extractSplicingSummary(aSwitchList, asFractionTotal = FALSE)
```

## Significance Thresholds

| Parameter | Default | Notes |
|-----------|---------|-------|
| isoform_switch_q_value | < 0.05 | Switch significance |
| dIF (delta isoform fraction) | > 0.1 | Minimum biological effect |
| Consequence q-value | < 0.05 | Significance per consequence type |
| Gene-level OFDR (stageR) | < 0.05 | Gene-level screening FDR |
| satuRn alpha | 0.05 | Empirical-Bayes alpha |
| swish qvalue | < 0.05 | Local FDR from qvalue package |

## Common Errors

Messages quoted from runs (IsoformSwitchAnalyzeR 2.6.0, DRIMSeq 1.34.0, fishpond 2.12.0); "(sic)" marks the packages' own spelling.

| Error | Cause | Solution |
|-------|-------|----------|
| `importRdata: The annotation and quantification ... seems to be different (Jaccard similarity < 0.925)` | Isoform IDs in the quantification do not match the GTF/FASTA. Version-suffixed IDs (`ENST...1`) against an unversioned GTF give "Only 0 overlap" | Pass `ignoreAfterPeriod = TRUE` to **both** `importIsoformExpression` and `importRdata` (verified: imported after that); otherwise rebuild the index from the same annotation |
| `In the designMatrix the following column(s): batch Contain constant information` | A covariate column with one value | Drop it |
| `The supplied design matrix will result in a model matrix that is not full rank` | Covariate identical to condition | Drop the covariate or re-design |
| `A statistical test cannot be performed without replicates` | One sample in a condition | Add replicates |
| `No genes were considered switching with the used cutoff values` | Test called with `reduceToSwitchingGenes = TRUE` and nothing passed (often the count route, see above) | Use `reduceToSwitchingGenes = FALSE` and check the count route |
| `Please run the 'addORFfromGTF()' ... to detect ORFs` | `extractSequence()` / `analyzePFAM()` before ORFs exist (GTF without CDS) | `analyzeORF()` first (only when there is no CDS) |
| `To test differences in <signal peptides / protein domains / IDR>, the result of the <SignalP / Pfam / NetSurfP2> analysis must be advailable (sic)` | Consequence type requested without its annotator imported | Drop the type or import the annotation |
| `The 'removeNoncodinORFs' argument must be supplied` | `analyzeCPC2()` called without it | Pass `FALSE` (see above) |
| `analyzePFAM: more columns than column names` | Raw `hmmscan --domtblout` file | Convert with `examples/hmmscan_to_pfamscan.py` |
| `The file(s) 'pathToSignalPresultFile' points to does not exist` (also `pathToIUPred2AresultFile`) | Annotator output missing or wrong path | Check the path; run the annotator |
| `dmFilter: !No genes left after filtering!` / `min_samps_gene_expr >= 0 && min_samps_gene_expr <= ncol(x@counts) is not TRUE` | Thresholds nothing passes / more samples than exist | Scale thresholds to the design |
| `swish: is.factor(condition) is not TRUE` | Character condition column | `factor()` it |
| `scaleInfReps: there are no inferential replicates in the assays of 'y'` | Salmon run without `--numGibbsSamples` | Re-run Salmon with `--numGibbsSamples 20` |
| `unable to find an inherited method for function 'samples' for signature 'object = "dmDSdata"'` | DEXSeq attached after DRIMSeq masks `samples()` | `DRIMSeq::samples(d)` |

## Common Pitfalls

- **Treating short-read-derived isoform calls as ground truth** -> Salmon EM is uncertain; use Gibbs samples + swish if quantification uncertainty matters.
- **Comparing across annotations** -> GENCODE basic vs comprehensive, RefSeq, Ensembl all have different transcript catalogs; switches "appear" or "disappear" with annotation choice. Document version.
- **Not running long-read where possible** -> Iso-Seq / ONT removes ambiguity for genes with many similar isoforms (TTN, MAPT, NEFM, DSCAM).
- **Reporting a "switch" without a sashimi plot** -> reviewers will demand it; do it upfront.
- **Forgetting stageR also corrects gene-level p when starting from DRIMSeq** -> DRIMSeq's gene-level adjusted p should be passed as `pScreen`, not raw transcript p-values.

## Related Skills

- differential-splicing - Event-level (rMATS, leafcutter, MAJIQ) complementary to DTU
- splicing-quantification - PSI is a 1D projection of DTU shifts
- splicing-qc - Verify upstream library, depth, alignment before DTU
- sashimi-plots - Required visualization for switch validation and reporting
- splice-variant-prediction - Connects SpliceAI variant predictions to specific isoforms
- long-read-splicing - Full-isoform DTU bypasses transcript-quant uncertainty; preferred for many-isoform genes
- pathway-analysis/go-enrichment - Pathway enrichment of switching genes
- rna-quantification/alignment-free-quant - Salmon with `--numGibbsSamples` is upstream

## References

- Han et al 2026 *NAR Genom Bioinform* 8(3):lqag098, 10.1093/nargab/lqag098 (preprint *bioRxiv* 10.64898/2025.12.08.693027) - IsoformSwitchAnalyzeR v2
- Vitting-Seerup & Sandelin 2019 *Bioinformatics* 35:4469-4471 - IsoformSwitchAnalyzeR original
- Anders et al 2012 *Genome Res* - DEXSeq
- Nowicka & Robinson 2016 *F1000Research* - DRIMSeq
- Gilis et al 2022 *F1000Research* - satuRn
- Zhu et al 2019 *NAR* - swish / fishpond
- Van den Berge et al 2017 *Genome Biol* - stageR
- Love, Soneson, Patro 2018 *F1000Research* - Swimming downstream DTU workflow
- Maquat 2004 *Nat Rev Mol Cell Biol* - NMD review
- Lykke-Andersen & Jensen 2015 *Nat Rev Mol Cell Biol* - NMD update
- Lindeboom et al 2016 *Nat Genet* - NMD rules and escape from RNA-seq
- Lareau et al 2007 *Nature* - SR protein AS-NMD autoregulation
- Ni et al 2007 *Genes Dev* - ultraconserved-element AS-NMD in splicing regulators
- Cuccurese et al 2005 *NAR* 33:5965-5977 - ribosomal protein AS-NMD autoregulation
- Brown et al 2022 *Nature* - UNC13A cryptic exon (TDP-43)
- Melamed et al 2019 *Nat Neurosci*; Klim et al 2019 *Nat Neurosci* - STMN2 premature polyadenylation (TDP-43)
- Han et al 2020 *Sci Transl Med* - SCN1A poison exon ASO
