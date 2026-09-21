---
name: bio-outlier-splicing-detection
description: Detects aberrant splicing in single rare-disease patients vs a control panel for research use, using FRASER 2 (Bioconductor; beta-binomial model on the Intron Jaccard Index with a PCA or autoencoder fit, default delta cutoff 0.1, q hyperparameter), OUTRIDER (gene-level outlier expression via autoencoder denoising), LeafcutterMD (Dirichlet-multinomial outlier mode of LeafCutter for annotation-free junctions), and DROP (Snakemake pipeline integrating FRASER2 + OUTRIDER + monoallelic expression). The statistical model is fundamentally different from differential splicing - single-sample-vs-cohort outlier detection rather than two-group comparison. Widely used in rare-disease research programs (Solve-RD, NIH UDN). Use when applying RNA-seq to undiagnosed Mendelian disease, prioritising predicted splice variants for follow-up, or detecting cryptic splicing in disease tissue. Not a clinical report or variant classification.
tool_type: r
primary_tool: FRASER
license: MIT
---

## Version Compatibility

Checked on (2026-09-20): FRASER 2.2.0 + OUTRIDER 1.24.0 (R 4.4.3, Bioconductor 3.20) and FRASER 2.6.1 + OUTRIDER 1.28.1 + DROP 1.6.1 (R 4.5.3, Bioconductor 3.22); LeafcutterMD from leafcutter 0.2.9 with regtools 1.0.0.

Install: `BiocManager::install(c('FRASER', 'OUTRIDER'))`; leafcutter (LeafcutterMD) from GitHub (`davidaknowles/leafcutter`, not Bioconductor); `conda install -c bioconda regtools snakemake star samtools bcftools`; DROP as in the DROP section.

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Outlier Splicing Detection

For rare-disease RNA-seq, the question is not "what differs between groups?" but "what is aberrant in this single patient relative to a panel of unaffected samples?". The statistical framework is **single-sample-vs-cohort outlier detection**, fundamentally different from two-group differential splicing.

## Research Use and Licences

- **Research use only.** Outliers are candidates for follow-up, not diagnoses. This Skill produces no clinical report, no ACMG/AMP evidence strength (PS3, PP3 or otherwise) and no pathogenicity call for an individual. Classification and reporting belong to a qualified clinical geneticist in an accredited laboratory with a validated assay.
- **Licences of the tools (this Skill's own `license: MIT` is what upstream declared).** FRASER 2.6.1, OUTRIDER 1.28.1 and DROP 1.6.1 are CC BY-NC 4.0 (academic, non-commercial; commercial use needs a licence from the authors; DROP prints the notice at `drop init`). The FRASER 2.2.0 and OUTRIDER 1.24.0 `LICENSE` files say MIT. Read the `LICENSE` file of the installed version.
- **Controls:** individual-level GTEx RNA-seq (BAMs) is controlled access (dbGaP); pooling it needs an approved application.

## Tool Taxonomy

| Tool | Statistic | Test target | Fails when |
|------|-----------|-------------|------------|
| FRASER 2 | Beta-binomial on Intron Jaccard Index; PCA fit (default) or autoencoder | Splicing outliers (per-sample, per-junction) | Cohort <20 samples; tissue mismatch |
| OUTRIDER | Autoencoder-denoised NB expression Z-score | Gene-level expression outliers (LoF, monoallelic) | Cohort <50 samples (see Cohort Size and Power) |
| LeafcutterMD | Dirichlet-multinomial outlier mode | Annotation-free intron usage (junction clusters only) | Few controls; intron retention / pseudoexon clusters give no p-value |
| DROP | Snakemake pipeline | All of above + monoallelic expression | Pipeline complexity for small projects |

Core reference: **FRASER 2** for splicing outliers, **OUTRIDER** for expression outliers, **DROP** to combine.

## Decision Tree by Scenario

| Scenario | Recommended approach |
|----------|----------------------|
| Single patient + panel of n>=50 controls | FRASER 2 (Intron Jaccard Index) |
| Single patient + small panel (n=20-50) | FRASER 2; choose q with `estimateBestQ`; tissue-matched, batch-matched controls |
| Patient + cohort <20 | Insufficient for outlier detection (missed at n=8-16); recruit more samples |
| Outlier expression suspected (loss of function, monoallelic) | OUTRIDER on a cohort of >=50 |
| Annotation-free outlier (cryptic exon, novel junction) | LeafcutterMD |
| Integrated pipeline (splicing + expression + MAE) | DROP |
| TDP-43 ALS post-mortem brain (cryptic exons) | FRASER 2; expect UNC13A, STMN2, ATG4B |
| SF3B1-mutant cancer sample | FRASER 2 with cohort-matched RNA-seq; expect cryptic 3'ss |
| Familial dysautonomia (ELP1) | FRASER 2 in fibroblast/iPSC; CNS tissue gives strongest signal |
| Stargardt deep-intronic ABCA4 | FRASER 2 in retina-relevant tissue |
| Solid tumor splicing biomarker | Differential splicing (n>=10 vs cohort) - see differential-splicing skill |
| RNA follow-up of a SpliceAI hit | FRASER 2 + cross-reference with predicted variant location (Variant + Outlier Integration) |

**Outlier vs differential:** outlier regime = single patient or small heterogeneous case series vs a control panel (single-sample p-value vs cohort distribution). Differential regime = two well-defined groups, n>=3 each (two-group test; see differential-splicing skill). With n>=10 patients sharing a phenotype prefer **differential** (more power).

## FRASER 2 Workflow

**Goal:** Detect aberrant splicing in patient samples vs cohort using the Intron Jaccard Index.

**Approach:** Count split reads per junction, compute the Intron Jaccard Index, fit a PCA (default) model to estimate expected values, then flag outliers by p-value and delta. Runnable version with argument handling: `examples/fraser2_rare_disease.R`.

```r
library(FRASER); library(BiocParallel)
bp <- if (.Platform$OS.type == 'windows') SerialParam() else MulticoreParam(8)

bam_files <- list.files('bams/', pattern = '\\.bam$', full.names = TRUE)
sample_table <- data.frame(
    sampleID = sub('\\.bam$', '', basename(bam_files)),
    bamFile = bam_files,
    pairedEnd = TRUE
)

# colData must be an S4Vectors DataFrame; a plain data.frame errors on FRASER 2.2.0 and 2.6.1
fds <- FraserDataSet(
    colData = S4Vectors::DataFrame(sample_table),
    workingDir = 'fraser_workdir',
    name = 'rare_disease_cohort'
)

fds <- countRNAData(fds, BPPARAM = bp)
fds <- calculatePSIValues(fds)
fds <- filterExpressionAndVariability(fds, minExpressionInOneSample = 20, minDeltaPsi = 0.0)

fitMetrics(fds) <- 'jaccard'
currentType(fds) <- 'jaccard'

# q: see Choosing q. FRASER 2.6.1: estimateBestQ(); FRASER 2.2.0 has only optimHyperParams()
fds <- estimateBestQ(fds, type = 'jaccard', plot = FALSE)
fds <- FRASER(fds, q = c(jaccard = bestQ(fds, 'jaccard')), implementation = 'PCA', BPPARAM = bp)

all_results <- as.data.frame(results(fds, psiType = 'jaccard', padjCutoff = 0.05, deltaPsiCutoff = 0.1))
patient_results <- all_results[all_results$sampleID == 'PATIENT_001', ]
patient_results <- patient_results[order(patient_results$padjust), ]
```

`filterExpressionAndVariability` is left at FRASER's variability defaults (`quantile = 0.75`, `quantileMinExpression = 10`).

**Differences from FRASER 1.x:** default metric is the single **Intron Jaccard Index** (1.x: psi5, psi3, theta); default `deltaPsiCutoff` **0.1** (1.x: 0.3); FRASER >=1.99.0 is FRASER 2. `implementation = 'PCA'` is the default fit; `'AE'` is an autoencoder (the old `correction=` argument is deprecated). Document the version and cutoff used.

**Reproducibility:** PCA fits are deterministic. AE fits differ run to run (measured on 2.6.1: 9,841 of 20,010 p-values differ >1e-6), and neither seed alone fixes it (`set.seed(1)` alone: 9,279 differ; `SerialParam(RNGseed = 1)` alone: 9,815 differ); `set.seed(1)` together with `BPPARAM = SerialParam(RNGseed = 1)` gives 0 differ. Use PCA unless you need AE.

## OUTRIDER for Gene-Level Outlier Expression

**Goal:** Detect genes with aberrantly high or low expression in patient samples.

**Approach:** Autoencoder denoising of expression matrix; outliers identified by Z-score and adjusted p-value. Needs **>=50 samples** (see Cohort Size and Power); at n<50 an empty result is expected, not an error.

```r
library(OUTRIDER); library(BiocParallel)

countTable <- read.table('counts.tsv', header=TRUE, row.names=1)
ods <- OutriderDataSet(countData = countTable)

ods <- filterExpression(ods, minCounts=TRUE, filterGenes=TRUE)
# OUTRIDER 1.24.0: estimateBestQ() returns the number q; 1.28.1: returns the object, q in metadata
q_best <- estimateBestQ(ods)
if (is(q_best, 'OutriderDataSet')) { ods <- q_best; q_best <- metadata(ods)[['optimalEncDim']] }
q_best <- max(q_best, 2)    # OUTRIDER() requires q > 1; the optimal-hard-threshold estimate can be 1
bp <- if (.Platform$OS.type == 'windows') SerialParam() else MulticoreParam(8)
ods <- OUTRIDER(ods, q = q_best, BPPARAM = bp)

res <- results(ods, padjCutoff = 0.05, zScoreCutoff = 0)
patient_outliers <- res[res$sampleID == 'PATIENT_001', ]
```

OUTRIDER (Brechtmann 2018 *Am J Hum Genet*) catches loss-of-function alleles producing transcript collapse, monoallelic effects, and tissue-inappropriate expression - complements splice outlier detection. `zScoreCutoff = 0` is OUTRIDER's and DROP's default; raise it (e.g. 2) only as a deliberate extra filter.

## LeafcutterMD for Annotation-Free Outlier Intron Usage

**Goal:** Detect outlier intron usage relative to a control panel without annotation dependence.

**Approach:** Run LeafcutterMD (LeafCutter's Dirichlet-multinomial outlier mode for Mendelian disease) against the control panel. `leafcutter_cluster_regtools.py` is in the leafcutter repository (`clustering/`), `leafcutterMD.R` in `scripts/`; neither is installed by `BiocManager`. BAMs must be indexed and carry XS strand tags (regtools writes strand `?` otherwise and the clustering step drops those junctions).

```bash
for bam in *.bam; do
    regtools junctions extract -a 8 -m 50 -s XS "$bam" -o "${bam%.bam}.junc"
done

ls *.junc > juncfiles.txt
python leafcutter_cluster_regtools.py -j juncfiles.txt -o leafcutter -m 50 -l 500000

leafcutterMD.R \
    --num_threads 4 \
    --output_prefix patient_outlier \
    leafcutter_perind_numers.counts.gz
```

Outputs: `patient_outlier_clusterPvals.txt` (cluster x sample), `_pVals.txt` (intron x sample; row names `chr:start:end:clu_N_strand` give the coordinates of a cluster), `_effSize.txt` (intron x sample). The p-values are **raw**; control the FDR yourself:

```r
p <- read.table('patient_outlier_clusterPvals.txt', header = TRUE, sep = '\t', check.names = FALSE)
long <- data.frame(cluster = rep(rownames(p), ncol(p)), sampleID = rep(colnames(p), each = nrow(p)),
                   p = unlist(p, use.names = FALSE))
long <- long[!is.na(long$p), ]
long$q <- p.adjust(long$p, method = 'BH')
hits <- long[long$q < 0.05, ]
```

LeafcutterMD (Jenkinson 2020 *Bioinformatics*) is junction-only: intron retention and pseudoexon events produced no cluster p-value in testing. Useful when FRASER's model fits poorly or when novel-junction sensitivity matters.

## DROP Pipeline (Integrated Workflow)

**Goal:** Run FRASER2 + OUTRIDER + monoallelic expression in one Snakemake pipeline.

**Approach:** DROP is distributed via **bioconda** (not PyPI). Install in a dedicated environment, initialise an empty project directory, fill in the sample annotation and `config.yaml`; the pipeline handles counting, autoencoding and reporting.

```bash
# Install via bioconda (DROP is not on PyPI)
mamba create -n drop_env -c conda-forge -c bioconda drop --override-channels
conda activate drop_env

mkdir my_diagnostic_run && cd my_diagnostic_run
drop init          # takes no project-name argument; run it in the empty directory

# Edit config.yaml (checked on DROP 1.6.1): fill sampleAnnotation, geneAnnotation, genome, root, htmlOutputPath;
#  aberrantSplicing: run: true   (implementation: PCA, FRASER_version: "FRASER2", padjCutoff 0.1, deltaPsiCutoff 0.1)
#  aberrantExpression: run: true (implementation: autoencoder)
#  mae: run: true                (needs RNA BAMs matched to a VCF)

snakemake aberrantSplicing --cores 16    # or aberrantExpression / mae; no --use-conda: DROP's rules have no conda: directives
```

`drop demo` downloads a public 10-sample demo project; its `snakemake aberrantSplicing` took about 85 minutes on 10 cores and reported no significant outlier (min padjust 1), as expected for 10 samples. The **MAE module** tests allelic imbalance on heterozygous SNPs called from RNA-seq with a negative-binomial test (not a z-score); useful for monoallelic LoF that splicing/expression outliers miss.

## Variant + Outlier Integration

**Goal:** Connect a candidate splice-altering DNA variant to an RNA-level outlier in the same sample, as a research prioritisation.

**Approach:** Cross-reference SpliceAI hits with FRASER2 outliers in the same sample. `delta_max` is the largest of the four SpliceAI delta scores (DS_AG, DS_AL, DS_DG, DS_DL).

```bash
bcftools query -f '%CHROM\t%POS\t%INFO/SpliceAI\n' spliceai_annotated.vcf > spliceai_raw.tsv
```

```r
library(dplyr)

raw <- read.table('spliceai_raw.tsv', sep = '\t', quote = '', col.names = c('chrom', 'pos', 'spliceai'),
                  stringsAsFactors = FALSE)
raw <- raw[raw$spliceai != '.', ]
# INFO = ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL, one comma-separated entry per ALT allele;
# unscored alleles are '.', which max() ignores
raw$delta_max <- vapply(strsplit(raw$spliceai, ',', fixed = TRUE), function(e)
    suppressWarnings(max(0, as.numeric(unlist(lapply(strsplit(e, '|', fixed = TRUE), `[`, 3:6))), na.rm = TRUE)),
    numeric(1))

fraser_hits <- read.table('fraser_results.tsv', header = TRUE, sep = '\t')
# 'chr21' vs '21' otherwise joins to zero rows without a warning
strip_chr <- function(x) sub('^chr', '', x)
variants <- mutate(raw, chrom = strip_chr(chrom))
fraser_hits <- mutate(fraser_hits, seqnames = strip_chr(seqnames))
stopifnot(any(variants$chrom %in% fraser_hits$seqnames))

confirmed <- variants %>%
    filter(delta_max >= 0.2) %>%
    inner_join(
        fraser_hits %>% filter(sampleID == 'PATIENT_001', padjust < 0.05),
        by = c('chrom' = 'seqnames'),
        relationship = 'many-to-many'
    ) %>%
    filter(abs(pos - start) < 1000 | abs(pos - end) < 1000)
```

A SpliceAI hit with a concordant FRASER2 outlier in the same sample is a candidate for functional follow-up. Whether it counts as functional evidence, and at what strength, is a clinical-laboratory decision (ClinGen SVI splicing recommendations, Walker 2023) outside this Skill. The 1 kb window is a loose heuristic.

## Cohort Size and Power

Measured 2026-09-20 on synthetic data (FRASER: planted events in a 30-sample cohort, 2x75 reads; OUTRIDER: `makeExampleOutriderDataSet(n = 2000 genes, freq = 1e-3)`, 3 seeds, padj<0.05).

**FRASER 2** (single patient + n-1 controls, PCA, q from `estimateBestQ`, FRASER 2.6.1): n=8, 12, 16: planted exon-skipping missed (0 calls); n=20, 25: flagged (3 junctions); n=30: all 4 planted events (skipping, cryptic donor, intron retention, pseudoexon) flagged on both FRASER 2.2.0 and 2.6.1; 0 false calls at every n. 4 real chrX samples: runs, 0 calls (n too small).

**OUTRIDER** (fraction of injected outliers detected; false calls 0-2 per 3 seeds):

| Samples | OUTRIDER 1.24.0 | OUTRIDER 1.28.1 |
|---------|-----------------|-----------------|
| 20 | 0/107 | 0/107 |
| 30 | 0/156 | 5/156 (3%) |
| 50 | 111/297 (37%) | 123/297 (41%) |
| 60 | 171/359 (48%) | 201/359 (56%) |
| 100 | 414/583 (71%) | 434/583 (74%) |
| 200 | 857/1203 (71%) | 921/1203 (77%) |

The q that `estimateBestQ` returns matters: on a 3000-gene synthetic matrix with a hidden batch (8 planted outliers) the OUTRIDER block above found 0/8 at n=30 on both versions and, at n=100, 4/8 with the 1.24.0 heuristic q=20 versus 2/8 with the 1.28.1 optimal-hard-threshold q (1, raised to 2; 1 false call). Fixed small q does not rescue n=30 (q=2 or 4: 0/156 on both versions) and q=2 at n=50 gave up to 39 false calls. Practical rule: OUTRIDER **>=50 samples** (useful power from ~100); FRASER **>=20** (more is better; the 4-event result above is n=30). Autoencoder covariance is learned across samples, so batch- and tissue-matched controls matter more than count.

| Cohort size | FRASER 2 | OUTRIDER |
|-------------|----------|----------|
| n < 20 | Missed at n=8-16; do not interpret a null | No power |
| n = 20-49 | Detected at n=20-30 in testing | Effectively none (0-3% at n=20-30) |
| n >= 50 | Recommended | ~40-55% at 50-60 |
| n >= 100 | Best calibration | ~70-75% |

GTEx-derived tissue-matched controls can supplement small in-house cohorts but introduce batch effects (and need dbGaP access); use only when in-house n < 30, keep them a minority, and document the pooling strategy. The tools do not correct for batch.

## Tissue Choice for Mendelian RNA-seq

| Tissue | Pros | Cons | Genes captured |
|--------|------|------|----------------|
| Whole blood (PAXgene) | Easy, standard | Globin contamination; many disease genes silent | ~70-80% of clinical genes |
| Fibroblast (skin biopsy) | Reasonable expression | Requires culture; senescence variability | ~75-85% |
| Muscle biopsy | Best for muscular dystrophy | Invasive | ~85-90% for muscle disorders |
| iPSC-derived neuron / cardiomyocyte | Disease-relevant tissue | Cost, variability | ~95% if differentiation works |
| Urine sediment | Non-invasive | Low yield | ~50-60% |

For UDN-style cases: blood first, then fibroblast if blood lacks expression of the candidate gene. **A negative blood RNA-seq does not rule out a candidate gene that is silent in blood** - verify gene expression with the GTEx tissue panel before committing to the tissue.

**Detecting a tissue/batch mismatch before interpreting.** Strict tissue matching is required, and a mismatched patient does not always look like "hundreds of outliers": on a synthetic sample with globally shifted splicing, FRASER PCA reported 0 calls (the fit absorbed it), FRASER AE q=10 reported 15, LeafcutterMD 9. So do not read a zero-call patient as tissue-matched. Compare calls per sample; a sample far above the cohort median is suspect:

```r
sort(table(factor(all_results$sampleID, levels = colnames(fds))), decreasing = TRUE)[1:5]
```

## Choosing q

`q` is the latent-space (bottleneck) dimension of the fit; too low leaves confounders in the outlier signal, too high absorbs real signal, and it must stay far below the sample count. Estimate it per cohort instead of hard-coding it:

```r
fds <- estimateBestQ(fds, type = 'jaccard', plot = FALSE)    # FRASER 2.6.1: optimal hard threshold, fast
bestQ(fds, 'jaccard')
# exhaustive injected-outlier grid (FRASER 2.2.0: optimHyperParams(fds, type = 'jaccard', q_param = ...)):
fds <- estimateBestQ(fds, type = 'jaccard', useOHT = FALSE, q_param = c(2, 5, 10, 15), plot = FALSE)
plotEncDimSearch(fds, type = 'jaccard', plotType = 'auc')    # plotType = 'loss' also works
```

`plotEncDimSearch` without `plotType` shows the OHT singular-value plot and returns NULL (with a warning) after `useOHT = FALSE`. Measured on the 30-sample synthetic cohort (PCA): q=1-3 flagged 4/4 planted events, q=5-10 3/4 (the 60%-usage cryptic donor was missed), q=15 1/4; `estimateBestQ` chose q=1 (FRASER 2.6.1, OHT) and q=2 (2.2.0 grid), both 4/4. A fixed q=10 (the old default) is wrong for cohorts this small. OHT prints "Optimal latent space dimension is smaller than 2 ... set to 2" when the cohort is too small or degenerate (seen at n=8-12). For OUTRIDER, `estimateBestQ(ods)` (see its block) is the fast single-q estimate; `findEncodingDim` / `estimateBestQ(ods, useOHT = FALSE)` grids many q values and is slow.

## Per-Tool Failure Modes

### FRASER 2: Tissue Mismatch

**Trigger:** Patient sample from a different tissue than the controls. **Mechanism:** the fit learns tissue-specific patterns; a mismatched patient becomes a global outlier or is absorbed. **Symptom:** many implausible calls (AE) or none (PCA); see Tissue Choice for the check. **Fix:** use only controls from the patient's tissue.

### OUTRIDER: Few Controls

**Trigger:** Cohort <50 samples. **Symptom:** no convergence warning; `Warning: No significant events` (measured n=8-30, 0 false calls). **Fix:** recruit >=50 matched controls; a z-score on log-CPM against the controls is a cruder fallback. Do not read an empty result as "no expression outliers".

### LeafcutterMD: Cluster Count Limits

**Trigger:** Very few clusters in the patient sample (low coverage or filtered out). **Mechanism:** the Dirichlet-multinomial fit is unstable with few observations. **Symptom:** inflated or deflated p-values; few significant calls. **Fix:** increase coverage; relax filtering (`-m 10` instead of 50); or switch to FRASER2.

### DROP: Snakemake Pipeline Failures

**Trigger:** Missing dependencies or incompatible R/Bioconductor versions; a killed run leaves a lock. **Symptom:** a step fails partway with a cryptic R error. **Fix:** install DROP in its own bioconda environment (see DROP Pipeline); after a killed run `snakemake --unlock`, then resume.

## Reconciliation: When Outlier Tools Disagree

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| FRASER2 sig, OUTRIDER not | Splicing change without expression collapse | Standard splicing outlier; report |
| OUTRIDER sig, FRASER2 not | Expression LoF without splicing change | Likely promoter / regulatory; not splicing |
| Both sig at same gene | LoF allele triggering NMD on splicing-disrupted transcript | Strong combined evidence; expect downstream |
| LeafcutterMD sig, FRASER2 not | Novel cryptic event not in annotation | High-priority novel finding; investigate |
| All tools null but biology suggests change | Underpowered cohort or wrong tissue | Verify gene expression in tissue; recruit larger cohort |

## Disease-Specific Expectations

| Condition | Expected outlier signature | Tissue |
|-----------|----------------------------|--------|
| ALS / FTD (TDP-43 loss) | Cryptic exons in UNC13A, STMN2, ATG4B | Post-mortem brain ONLY |
| SF3B1-mutant MDS / CLL / uveal melanoma | Aberrant 3'ss ~10-30nt upstream of canonical | Bone marrow / tumor tissue |
| Spinal muscular atrophy (untreated SMN2) | SMN exon 7 skipping | Fibroblast / iPSC-MN |
| Familial dysautonomia (ELP1 c.2204+6T>C) | ELP1 exon 20 skipping (>=99% in CNS, partial elsewhere) | iPSC-neuron > fibroblast > blood |
| Deep-intronic CFTR / USH2A / CEP290 | Pseudoexon inclusion | Cognate disease tissue (lung / retina) |
| Duchenne muscular dystrophy (DMD) | Out-of-frame exon skipping pattern | Muscle biopsy |
| Stargardt (ABCA4) deep-intronic | Pseudoexon in retina | Retinal organoid / iPSC-RPE |

For each, the gene must be expressed in the queried tissue. Verify with GTEx before assuming a negative result rules out the gene.

## Common Errors

Messages observed on FRASER 2.2.0/2.6.1, OUTRIDER 1.24.0/1.28.1, DROP 1.6.1, leafcutter 0.2.9.

| Error / message | Cause | Solution |
|-----------------|-------|----------|
| `FraserDataSet`: `assignment of an object of class "data.frame" is not valid for slot 'colData'` | Plain data.frame passed as `colData` | `S4Vectors::DataFrame(sample_table)` |
| `failed to open BamFile: failed to load BAM index` (FRASER); `Unable to open BAM/SAM index` (regtools) | BAM not indexed / index missing | `samtools index`; check `samtools idxstats` |
| `Optimal latent space dimension is smaller than 2 ... set to 2` (`estimateBestQ`) | Cohort too small or degenerate (seen at n=8-12) | Add matched controls; do not interpret a null |
| `OUTRIDER(ods, q = <OutriderDataSet>)`: `operations are possible only for numeric, logical or complex types` | OUTRIDER 1.28.1 `estimateBestQ` returns the object | Use the `is(q_best, 'OutriderDataSet')` step in the OUTRIDER block |
| `Please provide for q an integer greater than 1 ...` (`OUTRIDER()`) | OHT estimate of 1 on OUTRIDER 1.28.1 | `max(q_best, 2)` step in the OUTRIDER block |
| `Warning: No significant events` (OUTRIDER) | Cohort <50 (expected), or no outliers | See Cohort Size and Power |
| `estimateBestQ`: could not find function (FRASER 2.2.0) | Added after 2.2.0 | `optimHyperParams(fds, type = 'jaccard')` |
| LeafcutterMD: empty counts file, then `non-character argument` / `undefined columns selected` | Junction strand `?` (no XS tags), so clustering drops all junctions | Use BAMs with XS:A tags |
| `drop init my_diagnostic_run`: `Got unexpected extra argument` | `drop init` takes no argument | `mkdir` + `cd` + `drop init` |
| Variant integration returns 0 rows | VCF uses `1`, FRASER uses `chr1` (or vice versa); join is silent | `strip_chr()` step in the integration block |

## Common Pitfalls

- **Using bulk differential-splicing tools for n=1 vs cohort** - rMATS, regular leafcutter, SUPPA2 are not designed for this. Use FRASER2 / LeafcutterMD.
- **Ignoring tissue choice** - see Tissue Choice.
- **Forgetting batch effects** - combining in-house and external (GTEx) controls introduces sequencing batch confounding that the tools do not correct; see Cohort Size and Power.
- **Skipping the variant + outlier integration** - an RNA-only outlier without a DNA variant suggests cellular state or technical artifact; a DNA-only prediction without RNA confirmation is only a prediction.
- **Treating all FRASER2 outliers as pathogenic** - many are benign tissue-specific variation; filter against a disease-gene panel you supply and the patient's phenotype.
- **Reporting a clinical result** - this Skill's output is research-use only (Research Use and Licences).

## Quality Thresholds

| Metric | Recommendation | Source |
|--------|----------------|--------|
| Cohort size | FRASER n>=20 (min), >=50 (recommended); OUTRIDER n>=50 | Measured, Cohort Size and Power |
| FRASER 2 padj | < 0.05 | Standard |
| FRASER 2 delta Jaccard | >= 0.1 (default in v2) | Scheller 2023 *AJHG* |
| OUTRIDER padj | < 0.05 (zScoreCutoff 0 = tool default) | Brechtmann 2018 *AJHG* |
| LeafcutterMD | BH q < 0.05 over cluster x sample p-values | Measured |
| Sequencing depth | >=50M PE reads/sample | Standard for AS analysis |
| Tissue and batch match between patient and controls | Required | Critical for calibration |

## Related Skills

- splice-variant-prediction - SpliceAI / Pangolin for in-silico prediction; integration target
- differential-splicing - When testing multiple patients vs controls (>=10 vs cohort)
- splicing-qc - Library / depth / tissue prerequisites
- variant-calling/clinical-interpretation - ACMG/AMP framework (a clinical-laboratory step outside this Skill)
- workflows/clinical-trial-pipeline - Trial-grade RNA-seq pipelines

## References

- Mertes et al 2021 *Nat Commun* - FRASER 1.x
- Scheller et al 2023 *Am J Hum Genet* - FRASER 2.0 (Intron Jaccard Index)
- Brechtmann et al 2018 *Am J Hum Genet* - OUTRIDER
- Jenkinson et al 2020 *Bioinformatics* - LeafcutterMD
- Yepez et al 2021 *Nat Protocols* - DROP pipeline
- Cummings et al 2017 *Sci Transl Med* - RNA-seq for muscular dystrophy diagnostics
- Kremer et al 2017 *Nat Commun* - RNA-seq for mitochondrial disease
- Brown et al 2022 *Nature* - UNC13A cryptic exon (TDP-43 / ALS)
- Klim et al 2019 *Nat Neurosci* - STMN2 cryptic splicing (ALS)
- Darman et al 2015 *Cell Rep* - SF3B1 cryptic 3'ss
- Walker et al 2023 *Am J Hum Genet* - ClinGen SVI splicing recommendations
