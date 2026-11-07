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
| `To test differences in <signal peptides / protein domains / IDR>, the result of the <SignalP / Pfam / installed NetSurfP importer> analysis must be advailable (sic)` | Consequence type requested without its annotator imported | Drop the type or import the annotation |
| `The 'removeNoncodinORFs' argument must be supplied` | `analyzeCPC2()` called without it | Pass `FALSE` (see above) |
| `analyzePFAM: more columns than column names` | Raw `hmmscan --domtblout` file | Convert with `scripts/hmmscan_to_pfamscan.py` |
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
