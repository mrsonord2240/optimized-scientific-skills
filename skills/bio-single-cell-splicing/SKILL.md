---
name: bio-single-cell-splicing
description: Analyzes alternative splicing at single-cell resolution. The first decision is library chemistry — 10X 3' is fundamentally limited (RT primes from poly-A, R2 falls in 3' UTR, <0.1 junction read per cell per AS event). Plate-based full-length methods (Smart-seq3, FLASH-seq, VASA-seq, STORM-seq) and single-cell long-read (MAS-Iso-seq, scISOr-Seq2) are the chemistries that give per-cell isoform structure. Tools include MARVEL (R, Smart-seq integrated), BRIE2 (Bayesian PSI with regulatory features and ELBO_gain test), scQuint (junction-cluster, plate-based; not for 10X), SpliZ (annotation-free Z-score), Psix (graph-smoothness regulated AS), and Sierra (alternative polyadenylation, often confused with AS). Use when analyzing isoform usage in scRNA-seq, identifying cell-type-specific splicing, or determining whether scRNA-seq chemistry supports splicing analysis at all.
tool_type: mixed
primary_tool: MARVEL
license: MIT
author: GPTomics
---

## Version Compatibility

Reference examples checked on: MARVEL 2.0.5 (R 4.4, Seurat 5.5), BRIE2 2.3.0, scQuint 0.3.3, Psix 0.10.8, Sierra 0.99.27 (GitHub; there is no 1.0 release), regtools 1.0.0, leafcutter 0.2.9, anndata 0.12, scanpy 1.11, pandas 2.3. SpliZ (Nextflow pipeline) was not run: Nextflow is not installed on the checking machine.

Install (GitHub sources; BRIE2's PyPI sdist is broken):
```bash
pip install git+https://github.com/huangyh09/brie                 # BRIE2; `pip install brie` fails (sdist lacks requirements.txt)
pip install tensorflow tf_keras tensorflow-probability             # brie-quant needs all three, and TF_USE_LEGACY_KERAS=1
pip install --no-build-isolation git+https://github.com/lareaulab/psix   # numpy must already be installed
pip install git+https://github.com/songlab-cal/scquint             # pulls torch, pyro, snakemake
```
```r
remotes::install_github('wenweixiong/MARVEL')    # archived from CRAN 2025-10
remotes::install_github('VCCRI/Sierra')          # not on Bioconductor
```
SpliZ is run as a Nextflow pipeline (below); it has no `setup.py`, so `pip install git+...SpliZ` does not work.

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Single-Cell Splicing Analysis

The fundamental decision is **chemistry**, not tool. Most droplet 3' scRNA-seq cannot support transcriptome-wide splicing inference because reverse transcription primes from the poly(A) tail and most reads land in the 3' UTR — far from CDS-region splicing events. Plate-based full-length methods and single-cell long-read sequencing are the chemistries that give per-cell isoform structure across the gene body.

## The 10X 3' Problem (Quantified)

Three compounding mechanisms make 10X Chromium 3' (v3.1, GEM-X, v4) hostile to splicing:

1. **3' enrichment**: median fragment <1 kb from poly(A); >70% of unique reads fall within 3' UTR.
2. **Short R2 (~91 nt)**: each read straddles at most one junction; usually none, because R2 lands in 3' UTR.
3. **PCR concatemers and TSO artifacts**: pollute junction detection; UMI collapse is gene-level, not isoform-level.

**Quantitative estimate:** Only a small fraction of cassette exons sit close enough to the polyA site to be sampled by 3' chemistry (empirical estimates from APA/3'-end atlases — see Tian & Manley 2017 *Nat Rev Mol Cell Biol* for the 3' UTR isoform landscape). Effective junction read yield from 10X 3' is **<0.1 per cell per AS event** — vs the 5-10 needed for stable per-cell PSI. Most splicing analyses on 10X 3' data report artifacts.

**The 5' kit (10X 5' GEX) does not solve this** — it shifts capture from 3' UTR to 5' UTR / TSS-proximal regions. Marginal improvement; not a transcriptome-wide solution. Note that V(D)J recovery requires the **10X Chromium Single Cell Immune Profiling kit** (with TCR/BCR-specific enrichment), not 5' GEX alone — postdocs designing immune-repertoire experiments must use the dedicated V(D)J kit.

## Decision: Does the Chemistry Support Splicing Analysis?

| Chemistry | Splicing analysis viable? | Best alternative if no |
|-----------|----------------------------|--------------------------|
| 10X 3' (Chromium v3, GEM-X, v4, Flex) | No (transcriptome-wide); maybe near-3'-end events | Sierra for APA |
| 10X 5' GEX | Limited; near-5'-end events only | Sierra for alternative TSS; switch to MAS-Iso-seq |
| Smart-seq2 | Yes (full transcript) | MARVEL or BRIE2 |
| Smart-seq3 / Smart-seq3xpress | Yes + UMI molecule counting | MARVEL or BRIE2 |
| FLASH-seq | Yes (faster, cheaper Smart-seq3) | MARVEL or BRIE2 |
| VASA-seq | Yes + total RNA (incl. nascent, IR) | MARVEL with IR analysis |
| STORM-seq | Yes + total RNA + ribodepletion | MARVEL with IR analysis |
| MAS-Iso-seq + 10X 5' (PacBio Kinnex) | Yes — full isoforms per cell | FLAMES, scNanoGPS, IsoQuant, see long-read-splicing |
| scISOr-Seq2 (PacBio + 10X) | Yes — full isoforms with cell-typing | FLAMES, IsoQuant |
| ONT direct cDNA scRNA | Yes | FLAMES |
| ONT direct RNA scRNA | Yes + native modifications | FLAMES |

## Tool Selection Matrix

| Tool | Best for | Input | Strengths | Fails when |
|------|----------|-------|-----------|------------|
| MARVEL | Smart-seq plate-based and (v2+) 10X droplet unified workflow | Plate or droplet BAMs + Seurat | SE/A5SS/A3SS/MXE/RI/AFE/ALE; modality classification; native Seurat integration; v2 droplet support | R-only |
| BRIE2 | Plate-based with regulatory feature prior | Plate BAM + GFF3 events | Bayesian variational PSI + ELBO_gain test; principled uncertainty; CLI-driven (`brie-count`, `brie-quant`) | TensorFlow dependency; slow at scale |
| scQuint | Plate-based annotation-free junction-cluster quantification (validated on Smart-seq2) | STAR junctions across cells | Cluster-level junction usage; latent Dirichlet | Authors recommend AGAINST use on 10X 3'/5' data (3'-bias confounds); plate-based only |
| SpliZ | Annotation-free discovery of cell-state-associated splicing | STAR-aligned BAMs | Per-gene Z-score; no event database needed | Annotation-free = power tradeoff |
| Psix | Regulated AS along trajectories | PSI + mRNA-per-exon tables (events x cells) + low-dimensional cell space | Tests smoothness of PSI over the cell-state neighbourhood; robust to dropout | Needs a cell-state latent space (e.g. PCA) upstream |
| Sierra | APA in 10X 3' (NOT splicing) | 10X BAM + GTF | Peak-calling 3' ends; DEXSeq DTU on UTR isoforms | APA only; not for cassette exons |
| pseudobulk leafcutter / rMATS | Between-cell-type differential splicing | Aggregated BAMs | Bulk-level statistical power | Loses within-cluster heterogeneity |
| MAS-Iso-seq + FLAMES | Full-length single-cell isoforms | 10X 5' + PacBio Kinnex | Full isoforms per cell at scale | Cost; complex pipeline |

## Decision Tree by Goal

| Goal | Recommended approach |
|------|----------------------|
| "Will my 10X 3' data support splicing?" | No transcriptome-wide; consider Sierra for APA. Note: scQuint authors recommend against use on 10X data |
| Cassette exon analysis in cell types from Smart-seq2 | MARVEL with `ComputePSI` + `AssignModality` + `CompareValues` |
| Discover cell-state-associated splicing without an event database | SpliZ |
| Test regulated AS along developmental pseudotime | Psix |
| Per-cell PSI with uncertainty in low-coverage cells | BRIE2 |
| Differential splicing between two well-defined cell types | Pseudobulk leafcutter or rMATS on aggregated BAMs |
| APA (alternative polyadenylation, often confused with AS) | Sierra |
| Full-length single-cell isoforms at scale | MAS-Iso-seq + FLAMES (long-read) |
| Microexons (3-27 nt) | Long-read or aligner with low overhang (uLTRA, deSALT, or STAR `--alignSJoverhangMin 6 --alignSJDBoverhangMin 1` plus a strict mismatch filter; the usual 8/3 is too strict) |
| snRNA-seq (nuclei) — IR question | Library captures nuclear RNA enriched for incomplete splicing — interpret IR cautiously |

## MARVEL Plate-Based Workflow

**Goal:** Run a unified workflow from STAR junctions to cell-type-specific splicing calls.

**Approach:** Build a wide splice-junction count matrix (rows = junctions keyed by `coord.intron`, columns = cells), assemble per-event feature tables from rMATS, then construct MARVEL object with named slots (`SpliceJunction`, `SplicePheno`, `SpliceFeature`, `IntronCounts`, `GeneFeature`, `Exp`, `GTF`). Quantify PSI per event class, classify modality, test differential splicing.

Event tables (once per annotation): an rMATS run writes `fromGTF.<TYPE>.txt` (`rmats.py --task both`; one BAM per group is enough), and `Preprocess_rMATS` turns each into MARVEL's `tran_id` table. The GTF must be GENCODE-style (`gene` rows with a `gene_type` attribute; otherwise `gene_type` comes back NA) and the R object **must be named `gtf`**: `Preprocess_rMATS` reads that global and ignores its `GTF=` argument.

```r
library(MARVEL); library(data.table)

gtf <- fread('annotation.gtf', header=FALSE, sep='\t', quote='', data.table=FALSE)   # same GTF as the rMATS run
for (ev in c('SE', 'MXE', 'RI', 'A5SS', 'A3SS')) {
    tab <- Preprocess_rMATS(read.table(sprintf('rmats_out/fromGTF.%s.txt', ev), header=TRUE, sep='\t'),
                            GTF=gtf, EventType=ev)
    write.table(tab, sprintf('events_%s.txt', ev), sep='\t', quote=FALSE, row.names=FALSE)
}
```

```r
library(MARVEL); library(Seurat); library(data.table)

seurat_obj <- readRDS('cells.rds')   # meta.data: one row per cell, cell-type column 'cell.type'

# Build wide SJ matrix: first column 'coord.intron' (e.g. 'chr1:100007082:100022621'),
# subsequent columns are per-cell sample IDs with junction counts as values.
# This is constructed from STAR SJ.out.tab files (one per cell) merged on intron coord.
sj_files <- list.files('star_pass2/', pattern='SJ.out.tab$', full.names=TRUE)
sj_long <- rbindlist(lapply(sj_files, function(f) {
    d <- fread(f, sep='\t', header=FALSE,
               col.names=c('chr','start','end','strand','motif','annot','unique','multi','overhang'))
    d$coord.intron <- paste(d$chr, d$start, d$end, sep=':')
    d$sample <- gsub('_SJ.out.tab$', '', basename(f))
    d[, .(coord.intron, sample, unique)]
}))
sj <- dcast(sj_long, coord.intron ~ sample, value.var='unique', fill=0)

# SpliceFeature is a NAMED LIST keyed by event class (tables from Preprocess_rMATS above)
df.feature.list <- list(
    SE   = read.table('events_SE.txt',   header=TRUE, sep='\t'),
    A5SS = read.table('events_A5SS.txt', header=TRUE, sep='\t'),
    A3SS = read.table('events_A3SS.txt', header=TRUE, sep='\t'),
    MXE  = read.table('events_MXE.txt',  header=TRUE, sep='\t'),
    RI   = read.table('events_RI.txt',   header=TRUE, sep='\t')
)

# SplicePheno: per-cell metadata; sample.id column maps to SpliceJunction column names
df.pheno <- seurat_obj@meta.data
df.pheno$sample.id <- rownames(df.pheno)

marvel <- CreateMarvelObject(
    SpliceJunction = sj,
    SplicePheno    = df.pheno,
    SpliceFeature  = df.feature.list,
    GeneFeature    = read.table('gene_features.tsv', header=TRUE, sep='\t'),   # gene_id, gene_short_name, gene_type
    Exp            = read.table('tpm.tsv', header=TRUE, sep='\t'),   # non-log TPM; NOT row.names=1: first column must stay gene_id
    GTF            = fread('annotation.gtf', header=FALSE, sep='\t', quote='', data.table=FALSE)   # data.frame, no header (V1..V9)
)

marvel <- CheckAlignment(marvel, level='SJ')
marvel <- ComputePSI(marvel, CoverageThreshold=10, EventType='SE')   # one event class per call; PSI matrix in marvel$PSI$SE is 0-1
marvel <- CheckAlignment(marvel, level='splicing')
marvel <- CheckAlignment(marvel, level='gene')
marvel <- TransformExpValues(marvel, offset=1, transformation='log2', threshold.lower=1)

neurons <- df.pheno$sample.id[df.pheno$cell.type == 'neuron']
glia    <- df.pheno$sample.id[df.pheno$cell.type == 'glia']

# Modality is assigned per cell group (one call per group)
marvel <- AssignModality(marvel, sample.ids=neurons, min.cells=5, seed=1)
head(marvel$Modality$Results)

marvel <- CompareValues(
    marvel,
    cell.group.g1 = neurons, cell.group.g2 = glia,
    min.cells = 5, method = 'wilcox', method.adjust = 'fdr',
    level = 'splicing', event.type = 'SE'
)
res <- marvel$DE$PSI$Table[['wilcox']]   # p.val.adj = FDR; mean.g1, mean.g2, mean.diff = mean.g2 - mean.g1, all in PSI x 100
head(res[order(res$p.val.adj), c('tran_id', 'gene_short_name', 'mean.g1', 'mean.g2', 'mean.diff', 'p.val', 'p.val.adj')])
```

MARVEL 2.0.5 quirks (each reproduced):
- `RI` events need an `IntronCounts` matrix (`coord.intron` x cells) passed to `CreateMarvelObject`; the block above computes SE only.
- `ComputePSI(EventType='SE')` and `Preprocess_rMATS(EventType='SE')` fail with `non-character argument` / `replacement has 1 row, data has 0` when the SE table has no minus-strand event, and `ComputePSI` fails the same way when `coord.intron` chromosome names differ from `tran_id` (`1` vs `chr1`).
- `CoverageThreshold` above the read depth returns an all-NA PSI matrix.
- `method='dts'` needs the `twosamples` package, which MARVEL does not install.

For 10X droplet data, MARVEL v2+ provides `CreateMarvelObject.10x()` and `AnnotateSJ.10x()` constructors. Verify the exact API via `?CreateMarvelObject.10x` in installed MARVEL.

MARVEL classifies events into modalities (Song 2017 *Mol Cell*): included (PSI~1), excluded (PSI~0), bimodal (mixture at 0/1), middle (peaked ~0.5), multimodal. Bimodality usually reflects mixed cell states or stochastic monoallelic-like bursting. Mid-modality (peaked at 0.5) can be technical (mixed cells in a droplet) — confirm with full-length data.

## BRIE2 Bayesian PSI

**Goal:** Estimate per-cell PSI with informative regulatory-feature prior; test cell-state association via likelihood-ratio testing on covariate effects.

**Approach:** BRIE2 is a CLI-driven workflow (`brie-count` for read counting, `brie-quant` for variational inference + LRT). Get a GFF3 of splicing events, count cell-barcoded junction reads, then fit the model with covariate testing.

Events: use BRIE's precomputed exon-skipping GFF3 (human GENCODE v25, mouse GENCODE vM12; `SE.gold.gff3` strictest, `SE.most.gff3` more events) from `https://sourceforge.net/projects/brie-rna/files/annotation/`; `examples/sc_splicing_brie2.py` has `fetch_splicing_events()`. Align reads to the same or a close genome version. `briekit-event` does not work (its console script crashes with `ModuleNotFoundError: parseTables`).

```bash
# 1. Count splicing events per cell. sample_list.tsv: BAM path <tab> cell id; BAMs sorted AND indexed
#    (droplet data: -s possorted.bam -b barcodes.tsv.gz instead of -S; tags default to --cellTAG CB --UMItag UR)
brie-count \
    -a splicing_events.gff3 \
    -S sample_list.tsv \
    -o brie_counts/ \
    -p 16

# 2. Fit BRIE2 with LRT against the cell covariates (cell_metadata.tsv: cell id + numeric feature columns)
brie-quant \
    -i brie_counts/brie_count.h5ad \
    -c cell_metadata.tsv \
    -o brie_quant.h5ad \
    --interceptMode gene \
    --LRTindex All \
    --testBase null \
    --MCsize 3 \
    -p 16
```

`--interceptMode gene` fits a gene-specific intercept (recommended); `--LRTindex All` tests all covariates; `--testBase null` uses the null model as the LRT reference. `brie-quant` has no `--seed`: a repeat run gave the same significant set (33 = 33 events at fdr < 0.05) but per-cell Psi moved by up to 0.033. Verify exact flag set via `brie-quant -h` in installed BRIE2.

**Gene filter (silent).** `brie-quant` drops events below `--minCount 50 --minUniqCount 10 --minCell 30 --minMIF 0.001`; only its log says so (`Filtered out 31 genes ...`: 31 of 50 events on 130 real Smart-seq2 cells). Lower the three counts (`--minCount 10 --minUniqCount 3 --minCell 10` kept 21 of 50) when the output has far fewer events than the input.

```python
import scanpy as sc
import pandas as pd

a = sc.read_h5ad('brie_quant.h5ad')    # only events that passed the gene filter
# a.varm['ELBO_gain' | 'pval' | 'fdr' | 'cell_coeff'] are events x tested features, in cell_metadata.tsv column order (a.uns['Xc_ids'])
lrt = pd.DataFrame({'ELBO_gain': a.varm['ELBO_gain'][:, 0], 'pval': a.varm['pval'][:, 0],
                    'fdr': a.varm['fdr'][:, 0], 'cell_coeff': a.varm['cell_coeff'][:, 0]}, index=a.var_names)
hits = lrt[lrt.fdr < 0.05].sort_values('fdr')     # the sign of cell_coeff gives the direction
psi = a.layers['Psi']                              # cells x events, shrunken per-cell PSI; layers 'Psi_95CI' (interval width) and 'Z_std' (SD on the logit scale) give its uncertainty
```

BRIE2 (Huang & Sanguinetti 2021 *Genome Biol*) uses a sequence-derived feature prior (exon length, GC content, splice site strength, motif counts) to regularize PSI estimates in low-coverage cells. The LRT-based covariate test answers "is this event associated with cell state?" without requiring per-cell PSI accuracy. **Threshold `fdr`, not `ELBO_gain`:** on 100 pure-null events 12 had raw p < 0.05 and 5 had `ELBO_gain` > 3, but none had fdr < 0.05; on a planted set (20 x dPSI 0.7, 10 x 0.3, 70 null; 80 cells) fdr < 0.05 called 20/20, 10/10 and 3/70 null events, all with the planted sign.

## scQuint Junction-Cluster Differential Splicing

**Goal:** Test differential intron usage between two cell groups without an event database (plate-based data; not 10X).

**Approach:** Load STARsolo junction output (`--soloFeatures Gene SJ`; Smart-seq mode also gives it), annotate introns with a GTF, group introns that share a 3' splice site, then test each intron group.

```python
import numpy as np
from scquint.data import load_adata_from_starsolo, add_gene_annotation, group_introns
from scquint.differential_splicing import run_differential_splicing

adata = load_adata_from_starsolo('Solo.out/SJ/raw')        # dir with matrix.mtx, barcodes.tsv, features.tsv (SJ.out.tab columns)
adata = add_gene_annotation(adata, 'annotation.gtf.gz')     # Ensembl-style contigs (1, 2, X): scQuint prepends 'chr' to match the junctions
adata = group_introns(adata, by='three_prime')
adata.obs['cell_type'] = cell_types                         # obs holds only the barcodes: add your labels, one per cell
cell_idx_a = np.where(adata.obs.cell_type == 'neuron')[0]
cell_idx_b = np.where(adata.obs.cell_type == 'glia')[0]

intron_groups, introns = run_differential_splicing(adata, cell_idx_a, cell_idx_b,
                                                   min_cells_per_intron_group=10, min_total_cells_per_intron=10)
hits = intron_groups[intron_groups.p_value_adj < 0.05]      # introns: psi_a, psi_b, delta_psi per intron
```

The default `min_cells_per_intron_group=30` / `min_total_cells_per_intron=30` kept 10 of 200 planted intron groups on 60 cells with 22% low-coverage cells; when every group is filtered, both returned tables are empty. Lowered to 10, the planted 200 events gave 40/40 and 20/20 planted events and 2/140 null events at `p_value_adj` < 0.05.

## SpliZ for Annotation-Free Discovery

**Goal:** Identify splicing-defined cell populations without an event database.

**Approach:** Compute per-gene splicing Z-score across cells; test for cell-state association via permutation.

```bash
# SpliZ is a Nextflow pipeline (not a standalone CLI). Configure inputs in a .config
# file (dataname, input_file, libraryType, grouping_level_1/2) - either SICILIAN
# output (SICILIAN=true) or BAMs via a samplesheet CSV + metadata + GTF (SICILIAN=false).
nextflow run salzmanlab/spliz -r main -latest -c spliz.config
```

Not run here (see Version Compatibility); the config keys above match the repository's `nextflow.config`.

SpliZ (Olivieri 2022 *Nat Methods*) is robust to dropout because it pools junction information across the gene; particularly useful for discovering splicing diversity in heterogeneous tumor samples.

## Psix for Regulated AS Along Trajectories

**Goal:** Detect AS that varies coherently with cell state along a developmental trajectory, robust to dropout.

**Approach:** Score whether observed PSI is smooth over cell neighbourhoods defined in a low-dimensional expression space. Psix builds its own neighbour metric from a latent space you supply; it does not use a `sc.pp.neighbors` graph or `adata.obsp`.

```python
import pandas as pd
import scanpy as sc
import psix

adata = sc.read_h5ad('cells.h5ad')                       # normalised expression with X_pca
latent = pd.DataFrame(adata.obsm['X_pca'][:, :20], index=adata.obs_names)   # cells x dims (or a TSV path)

# psi_matrix.tsv, mrna_matrix.tsv: events x cells, first column = event id, header = cell ids.
# PSI may contain NaN; mrna = estimated mRNA molecules captured per event per cell (from TPM, or UMI counts).
psix_obj = psix.Psix(psi_table='psi_matrix.tsv', mrna_table='mrna_matrix.tsv')
psix_obj.run_psix(latent=latent, n_jobs=4)               # defaults n_random_exons=2000, n_neighbors=100

regulated = psix_obj.psix_results.query('qvals < 0.05')   # columns: psix_score, pvals, qvals
```

To build the two tables from per-cell STAR `SJ.out.tab` files, `psix.Psix().junctions2psi(sj_dir=..., intron_file=<Psix cassette-exon annotation>, tpm_file=..., save_files_in='psix_out/')` writes `psi.tab.gz` and `mrna.tab.gz` (Psix README; not run here). On 300 planted trajectory cells x 150 exons (30 dynamic) `qvals < 0.05` called 30/30 dynamic and 1/120 static exons with `n_random_exons=300, n_neighbors=30` (4/120 with the defaults, which took about 4 min on 4 jobs); median `psix_score` was 1.85 for dynamic and -0.04 for static exons.

Psix (Buen Abad Najar 2022 *Genome Res* 32:1385) is the principled alternative to imputing PSI: do not impute (it obliterates heterogeneity); test for graph smoothness instead.

## Sierra for APA (Not Splicing)

**Goal:** Detect alternative polyadenylation in 10X 3' data — frequently confounded with AS.

**Approach:** Peak-call read pile-ups at 3' ends, then DEXSeq-style DTU on 3' UTR isoforms. `FindPeaks` needs a junction file.

```bash
# regtools 1.0.0: -s is required; 10X R2 reads are sense-strand, so -s RF. (-s XS writes strand '?' on BAMs without XS tags.)
regtools junctions extract -a 8 -m 50 -M 500000 -s RF -o junctions.bed possorted_genome_bam.bam
```

```r
library(Sierra)

FindPeaks(
    output.file = 'peaks.txt',
    gtf.file = 'annotation.gtf',
    bamfile = 'possorted_genome_bam.bam',
    junctions.file = 'junctions.bed'      # regtools BED or STAR SJ.out.tab
)

# CountPeaks writes a MEX directory and returns NULL; read it back
CountPeaks(
    peak.sites.file = 'peaks.txt',
    gtf.file = 'annotation.gtf',
    bamfile = 'possorted_genome_bam.bam',
    whitelist.file = 'barcodes.tsv',
    output.dir = 'peak_counts/'
)
counts <- ReadPeakCounts(data.dir = 'peak_counts/')     # peak x cell sparse matrix

# AnnotatePeaksFromGTF also writes a file (returns NULL)
AnnotatePeaksFromGTF(
    peak.sites.file = 'peaks.txt',
    gtf.file = 'annotation.gtf',
    output.file = 'peak_annotations.txt'
)
peak.annotations <- read.table('peak_annotations.txt', header = TRUE, sep = '\t',
                               row.names = 1, stringsAsFactors = FALSE)

# cell_identities: named vector, barcode -> population. Use NewPeakSCE: DUTest on a
# NewPeakSeurat object fails with SeuratObject >= 5 (GetAssayData `slot` is defunct).
# min.cells/min.peaks default to 10/200: a small peak set then drops every cell and DUTest errors (invalid 'row.names' length)
peaks.sce <- NewPeakSCE(peak.data = counts, annot.info = peak.annotations, cell.idents = cell_identities,
                        min.cells = 0, min.peaks = 0)

apa_results <- DUTest(peaks.sce, population.1 = 'ctrl', population.2 = 'trt')   # gene_name, padj, Log2_fold_change per peak
```

`DUTest` pools cells into `num.splits = 6` random pseudobulk profiles per population; pass `replicates.1` / `replicates.2` (lists of cell sets, e.g. per donor) when donors exist. On the bundled BAM `CountPeaks` UMI totals were within 1% of an independent pysam count of distinct (CB, UB) pairs per peak, a random two-way split returned no peaks, and thinning one peak in each of 6 genes to 20% in one population was recovered 6/6 (1 of 5 unchanged genes was also called).

If only 10X 3' data is available, this is often what is actually wanted. Distinct UTRs change miRNA targeting, RBP binding, and stability — biologically meaningful but not splicing.

## Pseudobulk for Statistical Power

**Goal:** Recover bulk-level statistical power for differential splicing between cell types.

**Approach:** Sum junction counts per **sample x cell type** (donor, batch or plate), not per cell type: pooling all cells of a type into one column is n = 1 per group, and a Fisher test on the pooled counts called 34 of 140 planted null events at FDR < 0.05. Require >= 3 replicates per group, then run leafcutter on the pseudobulk columns.

```python
import numpy as np
import pandas as pd

def pseudobulk_junctions(junction_counts, cell_metadata, groupby='cell_type',
                         sample_col='sample', min_replicates=3):
    """junction_counts: junctions x cells. cell_metadata: one row per cell, INDEXED BY CELL ID,
    with columns groupby and sample_col. Returns (counts, groups): junctions x (cell type, sample)
    pseudobulks and a leafcutter groups table."""
    missing = junction_counts.columns.difference(cell_metadata.index)
    if len(missing):
        raise ValueError(f'{len(missing)} cells not in cell_metadata.index (set index=cell id): {list(missing[:3])}')
    meta = cell_metadata.loc[junction_counts.columns, [groupby, sample_col]]
    if meta.isna().any().any():
        raise ValueError(f'NaN in {groupby!r} or {sample_col!r}: those cells would be dropped')
    pb = junction_counts.T.groupby([meta[groupby].to_numpy(), meta[sample_col].to_numpy()]).sum().T
    assert pb.to_numpy().sum() == junction_counts.to_numpy().sum()   # every read kept
    groups = pd.DataFrame({'sample': [f'{g}__{s}' for g, s in pb.columns],
                           'group': [g for g, _ in pb.columns]})
    n_rep = groups.group.value_counts()
    if (n_rep < min_replicates).any():
        raise ValueError(f'need >= {min_replicates} samples per {groupby}; got {n_rep.to_dict()}')
    pb.columns = groups['sample']
    return pb, groups
```

Cells are joined to the metadata by id, and a metadata table with a default `RangeIndex` raises (a positional join would return all-zero columns with no error).

Feed leafcutter with junction rows already in cluster form (`chr:start:end:clu_N_strand`, the `*_perind_numers.counts.gz` of `leafcutter_cluster_regtools.py`, read with `pd.read_csv(..., sep=' ', index_col=0)`):

```python
pb, groups = pseudobulk_junctions(junction_counts, cell_metadata)
pb.to_csv('pb_counts.txt.gz', sep=' ')
groups.to_csv('pb_groups.txt', sep='\t', header=False, index=False)   # two groups only: subset the two cell types for each pair
```
```bash
leafcutter_ds.R -i 5 -g 3 -c 20 -o pb_ds pb_counts.txt.gz pb_groups.txt
```

leafcutter defaults are `-i 5` (min samples using an intron), `-g 3` (min samples per group) and `-c 20` (min reads per sample); with fewer replicates it tests nothing unless lowered, and `-g 1 -i 1` on one pseudobulk per type reports p-values with no replicate support. On 6 v 6 planted pseudobulks (5 cells each) it called 40/40 and 20/20 planted events and 1/140 null events at FDR < 0.05 (16/140 at raw p < 0.05: cells are overdispersed, so keep the FDR cut). rMATS takes BAMs, not count tables: merge each sample x cell type's cell BAMs with `samtools merge` first (not run here).

Use pseudobulk for differential splicing **between** well-defined cell types; use per-cell methods for **within-population heterogeneity** (graded splicing along pseudotime, bimodal cell-state mixtures).

## Single-Cell Long-Read = Future of Single-Cell Splicing

In 2024-2026, full-length single-cell long-read sequencing has become practical and is the recommended chemistry for splicing-focused single-cell experiments:

- **MAS-Iso-seq / PacBio Kinnex**: concatenated full-length cDNA arrays, ~16x throughput vs plain Iso-Seq, compatible with 10X 5' libraries (Al'Khafaji 2024 *Nat Biotech*)
- **scISOr-Seq2**: hybrid 10X + PacBio for cell typing + isoform structure (Joglekar et al 2024 *Nat Neurosci* 27:1051-1063, single-cell long-read brain isoform mapping)
- **ONT direct cDNA + 10X**: lower cost, similar information content
- **FLAMES**: barcode demultiplexing + isoform quantification + SNV calling for ONT scRNA (Tian 2021 *Genome Biol* 22:310)

For splicing-specific full-length single-cell analysis, see `long-read-splicing` skill (align with minimap2 `splice:hq`, quantify with Bambu or IsoQuant).

## Per-Tool Failure Modes

### MARVEL: SpliceJunction Matrix Format

**Trigger:** Building the SpliceJunction matrix from STAR SJ.out.tab incorrectly (e.g. long-format instead of wide).

**Mechanism:** MARVEL plate-based `CreateMarvelObject(SpliceJunction = ...)` expects a **wide matrix** with first column `coord.intron` (formatted `chr:start:end`) and subsequent columns being per-cell sample IDs with integer junction counts. Long-format data.frames or missing `coord.intron` column cause runtime errors.

**Symptom:** `$ operator is invalid for atomic vectors` (long format) or `undefined columns selected` (first column not named `coord.intron`); or an empty/all-NA PSI table despite junction reads (chromosome naming differs from the event `tran_id`, or `CoverageThreshold` above the depth).

**Fix:** Verify wide-matrix structure; ensure SJ.out.tabs are merged on the `chr:start:end` key with cells as columns. Use `data.table::dcast` for the long->wide reshape.

### BRIE2: TensorFlow Memory

**Trigger:** Large cohort (>10k cells) with deep coverage.

**Mechanism:** Variational inference loads full count matrix; TensorFlow allocates GPU memory aggressively.

**Symptom:** OOM kills; training stalls.

**Fix:** Reduce `--batchSize` from default (500000) to 100000 or 50000; train per-chromosome batch; use CPU mode for very small cohorts. Note flag is camelCase `--batchSize`, not `--batch_size`.

### scQuint: 3' Data Sparsity

**Trigger:** Running scQuint on 10X 3' v3 data hoping for splicing signal.

**Mechanism:** scQuint filters introns and intron groups by the number of cells with reads (defaults 30 per group) before testing; 10X 3' yields too few junction reads to pass.

**Symptom:** Most or all intron groups are filtered (10 of 200 survived on 60 sparse plate cells); when none survive, both result tables are empty.

**Fix:** Pivot to APA analysis with Sierra; or upgrade chemistry to MAS-Iso-seq. On plate data with few cells lower `min_cells_per_intron_group` and `min_total_cells_per_intron`.

### Psix: Wrong Inputs

**Trigger:** Building Psix from an AnnData, a `sc.pp.neighbors` graph or a `psi_matrix_path`.

**Mechanism:** Psix takes PSI and mRNA tables (events x cells) plus a latent space (cells x dims) and computes its own cell metric; it never reads `adata.obsp['connectivities']`.

**Symptom:** `TypeError: Psix.__init__() got an unexpected keyword argument 'psi_matrix_path'`; `UndefinedVariableError: name 'pvalue'` when querying results (columns are `psix_score`, `pvals`, `qvals`).

**Fix:** Use the call in the Psix section.

### Sierra: Annotation Gaps

**Trigger:** GTF missing 3'UTR annotations.

**Mechanism:** Sierra peak-calls within annotated 3'UTRs; missing annotations mean missed peaks.

**Symptom:** Few peaks detected; gene-level coverage but no APA calls.

**Fix:** Use comprehensive GENCODE annotation; or run de-novo peak calling first.

## Reconciliation: When Single-Cell Tools Disagree

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| MARVEL sig, BRIE2 not | Per-cell PSI noise (BRIE2 conservative); MARVEL pseudobulk-like | Trust MARVEL for cell-type comparisons; BRIE2 for within-cluster |
| BRIE2 sig, MARVEL not | Cell-state effect smoother than cell-type boundary | Test along trajectory with Psix |
| SpliZ sig, MARVEL not | Annotation-free SpliZ catches novel events | Investigate junction structure manually |
| Sierra sig, MARVEL not | Sierra is APA, MARVEL is splicing — different biology | Distinguish in interpretation |
| Pseudobulk sig, per-cell not | Power issue; effect averaged out per-cell | Report at cluster level, not per-cell |

## Quantitative Concepts Unique to Single-Cell

**Per-cell PSI vs pseudobulk PSI:**
- Per-cell PSI: meaningful only above the junction-read threshold in Quality Thresholds (plate-based or long-read).
- Pseudobulk PSI: aggregate, recovers bulk-level statistical power, discards within-cluster heterogeneity.

**Modality detection in PSI distributions** (Song 2017 *Mol Cell*):
| Modality | PSI distribution | Biology |
|----------|------------------|---------|
| Included | Peaked at 1 | Constitutive inclusion |
| Excluded | Peaked at 0 | Constitutive skipping |
| Bimodal | Mixture at 0 and 1 | Mixed cell states or monoallelic-like bursting |
| Middle | Peaked ~0.5 | Often technical (well-contamination, doublets, or low-coverage shrinkage to prior); confirm with full-length |
| Multimodal | Multiple peaks | Complex regulation; deserves follow-up |

**Beta-binomial vs binomial models:** with sparse counts, binomial PSI is overdispersed. Beta-binomial models (BRIE2; leafcutter2 as Dirichlet-multinomial cluster-level) handle this. For very sparse droplet data, even beta-binomial fits poorly per cell — collapse to pseudobulk.

**Imputation pitfalls:** naive imputation (MAGIC, scImpute, ALRA) of expression matrices is **not** appropriate for PSI: imputing missing junction counts averages over neighboring cells and obliterates the very heterogeneity under study. Psix's approach — testing smoothness of observed PSI over cell neighbourhoods — is the principled alternative.

## Cell-Type-Specific Splicing Biology

| System | Event | Regulator |
|--------|-------|-----------|
| Neural microexons | 3-27 nt exons enriched in brain | SRRM4/nSR100 (Irimia 2014 *Cell*); SRRM3 in retina/photoreceptors (Ciampi 2022 *PNAS*) |
| Neural differentiation | PTBP1 -> PTBP2 switch | miR-124 represses PTBP1; derepresses neural exons (Boutz 2007 *Genes Dev*) |
| T-cell activation | CD45 RA -> RO | hnRNP-L, ESRP-mediated |
| Erythropoiesis | EPB41 exon 16 | Splicing factor switching during maturation |
| Cardiac development | TTN N2BA -> N2B | MBNL1/CELF1 antagonism |
| EMT | FGFR2 IIIb -> IIIc, ENAH exon 11a | ESRP1/2 loss in mesenchymal state (Warzecha 2009 *Mol Cell*) |
| Activated T cell | CD45 isoform shift | Multiple SR/hnRNP regulators |

## Quality Thresholds

| Metric | Recommendation |
|--------|----------------|
| Cells per event with reads | >=50 (per-cell PSI) |
| Junction reads per event per cell | >=5-10 for stable per-cell PSI (MARVEL `CoverageThreshold=10`); <=1 = unreliable |
| PSI variance for cell-type call | <0.1 within cluster, >0.2 between clusters |
| Library | full-length plate or long-read for transcriptome-wide; 3' for APA only |
| Doublet filtering | Required before splicing analysis (DoubletFinder, Scrublet) |
| Cells per cluster (pseudobulk) | >=100 ideal; >=50 minimum; >=3 replicate samples per group |
| nuclear vs whole-cell | snRNA-seq enriches IR; treat with caution |

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| MARVEL `ComputePSI` / `Preprocess_rMATS`: `non-character argument`, `replacement has 1 row, data has 0` | An event table with no minus-strand event, or `coord.intron` chromosome names differ from `tran_id` | Same naming in both; see MARVEL quirks |
| MARVEL `CheckAlignment(level='gene')`: `undefined columns selected` | `Exp` read with `row.names=1`, so `gene_id` is gone | Read `tpm.tsv` without `row.names` |
| MARVEL `CompareValues`: `there is no package called 'twosamples'` | `method='dts'` | Install `twosamples` or use `method='wilcox'` |
| `brie-count`: `Could not retrieve index file`, then `IndexError: list index out of range` | BAMs have no `.bai` | `samtools index` every BAM |
| Sierra `FindPeaks`: `unused argument (bam.file = ...)` | The argument is `bamfile`, and `junctions.file` is required | Use the call in the Sierra section |
| Sierra `DUTest`: ``The `slot` argument of `GetAssayData()` was deprecated ... and is now defunct`` | `NewPeakSeurat` object with SeuratObject >= 5 | Build the object with `NewPeakSCE` |

## Related Skills

- single-cell/preprocessing - QC and normalization (must run before splicing)
- single-cell/clustering - Cell type annotation prerequisite
- single-cell/doublet-detection - Doublet filtering critical for splicing
- single-cell/data-io - h5ad / Seurat I/O
- splicing-quantification - Bulk RNA-seq comparison context
- long-read-splicing - Full-isoform analysis from MAS-Iso-seq, scISOr-Seq2; future of single-cell splicing

## References

- Huang & Sanguinetti 2021 *Genome Biol* - BRIE2
- Wen et al 2023 *Nucleic Acids Research* 51:e29 - MARVEL
- Benegas, Fischer & Song 2022 *eLife* - scQuint (annotation-free single-cell splicing analysis, validated on Smart-seq2)
- Olivieri et al 2022 *Nat Methods* - SpliZ
- Buen Abad Najar et al 2022 *Genome Research* 32:1385 - Psix
- Patrick et al 2020 *Genome Biol* - Sierra
- Song et al 2017 *Mol Cell* - splicing modality classification
- Picelli et al 2014 *Nat Protoc* - Smart-seq2
- Hagemann-Jensen et al 2020 *Nat Biotech* - Smart-seq3
- Hagemann-Jensen et al 2022 *Nat Biotech* - Smart-seq3xpress
- Hahaut et al 2022 *Nat Biotech* - FLASH-seq
- Salmen et al 2022 *Nat Biotech* - VASA-seq
- Johnson et al 2022 *bioRxiv* 10.1101/2022.03.14.484332 - STORM-seq (preprint)
- Al'Khafaji et al 2024 *Nat Biotech* - MAS-Iso-seq / Kinnex
- Tian et al 2021 *Genome Biology* 22:310 - FLAMES
- Joglekar et al 2024 *Nat Neurosci* 27:1051-1063 - scISOr-Seq2 single-cell brain isoform mapping
- Irimia et al 2014 *Cell* - neural microexons / SRRM4
- Ciampi et al 2022 *PNAS* 119:e2117090119 - SRRM3-dependent photoreceptor microexons
- Boutz et al 2007 *Genes Dev* - PTBP1/PTBP2 neural switch
- Tian & Manley 2017 *Nat Rev Mol Cell Biol* - alternative polyadenylation and 3' UTR isoforms
