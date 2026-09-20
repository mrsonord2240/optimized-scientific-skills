---
name: bio-differential-splicing
description: Detects differential alternative splicing between conditions using rMATS-turbo (binomial LRT on junction counts), leafcutter (Dirichlet-multinomial GLM on intron clusters), MAJIQ V3 deltapsi/HET (Bayesian posterior on LSVs), SUPPA2 (empirical-null on TPM-derived PSI), or Shiba (junction-imbalance-corrected). Reports FDR-corrected significance and delta PSI effect sizes. Tools differ in statistical model, annotation dependence, calibration regime, and replicate-count requirements. Use when comparing splicing patterns between treatment groups, tissues, or disease states.
tool_type: mixed
primary_tool: rMATS-turbo
license: MIT
---

## Version Compatibility

Reference examples tested with: rMATS-turbo 4.3+, SUPPA2 2.4+, leafcutter 0.2.9+, MAJIQ 3.0+, Shiba 0.5+, STAR 2.7.11+, regtools 1.0+, pandas 2.2+, R 4.4+

Commands checked 2026-09-20 on rMATS-turbo 4.4.0, PAIRADISE 1.0, SUPPA2 2.4 (statsmodels 0.14.6), leafcutter 0.2.9, Shiba 0.8.2, regtools 1.0.0, samtools 1.24, statsmodels 0.15.0 (confounder snippet). **MAJIQ / VOILA are licence-gated and were not run**; every MAJIQ statement below is from its public documentation (MAJIQ 3.0.11.dev7) and unverified.

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Differential Splicing

Detect splicing changes between conditions. Tool choice is a decision about **statistical model**, **annotation dependence**, and **calibration regime** under the specific experimental design — not a preference. Wrong tool for the design produces uncalibrated FDR or systematic effect-size bias.

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

## Install Notes

- rMATS-turbo, regtools, SUPPA2, Shiba: `conda install -c bioconda rmats regtools suppa shiba` (rMATS is not on PyPI; SUPPA2 needs statsmodels <=0.14). `--paired-stats` also needs PAIRADISE (rMATS section).
- leafcutter: not on bioconda; GitHub (`devtools::install_github('davidaknowles/leafcutter/leafcutter')`). 0.2.9 needed patches to build against current rstan (Stan array syntax, TBB) and `doMC` is Unix-only, so use Linux/WSL. `leafcutter_ds.R` and the clustering scripts come from a clone of the repo.
- DESeq2/DEXSeq/limma (custom designs): Bioconductor.
- MAJIQ V3 / VOILA: academic-licence download from majiq.biociphers.org.

## rMATS-turbo Differential Analysis

**Goal:** Detect statistically significant differential splicing between two groups from BAMs.

**Approach:** Run rMATS-turbo without `--statoff`, then filter by FDR + ΔPSI + per-replicate coverage.

```bash
# -t paired = paired-end reads (use -t single for single-end); --readLength = the read length
rmats.py \
    --b1 condition1_bams.txt \
    --b2 condition2_bams.txt \
    --gtf annotation.gtf \
    -t paired \
    --readLength 150 \
    --variable-read-length \
    --libType fr-firststrand \
    --nthread 8 \
    --od rmats_output \
    --tmp rmats_tmp \
    --novelSS \
    --cstat 0.05
```

`--cstat 0.05` tests `|ΔPSI| > 0.05`; raise to 0.10 for stricter discovery. `--novelSS` enables novel-junction discovery (recommended with STAR 2-pass).

**Match the command to the data.** rMATS exits 0 and writes header-only tables when it cannot use the reads:
- `-t paired` on single-end BAMs gave 0 events; use `-t single` for single-end data.
- `--readLength` = the read length. On 75 nt reads, `--readLength 150` without `--variable-read-length` gave 0 events; with `--variable-read-length` it only sets IncFormLen/SkipFormLen. `examples/diff_splicing_rmats.sh` reads type and length from the first BAM.
- After every run assert that `SE.MATS.JC.txt` has rows and an `FDR` column before filtering.
- Strandedness: `--libType` must match the library (RSeQC `infer_experiment.py`); `fr-firststrand` on unstranded data cut SE events 226 -> 121.

**`--paired-stats`** (paired designs: entry *i* of the comma-separated `--b1` list pairs with entry *i* of `--b2`) needs the R package **PAIRADISE** in the R that rMATS finds. Without it rMATS logs `no package called PAIRADISE`, exits 0 and writes header-only tables with no FDR column, so run the table check above. Install (public, no login): `git clone https://github.com/Xinglab/PAIRADISE`, install R packages `nloptr doParallel foreach iterators`, then `R CMD INSTALL PAIRADISE/pairadise/src/pairadise_model`.

Progress is written to `<od>/tmp/JC_SE/pairadise_status.txt`. If it stalls at a few percent with idle R workers (seen on WSL2: PAIRADISE's default socket cluster hung in 3 of 6 runs), force forked workers: put `parallel:::setDefaultClusterOptions(type = "FORK")` in a file and run rMATS with `R_PROFILE_USER=<file>`. Checked on 5 vs 5 pairs (120 SE events): 31 s, 22/24 strong planted events found with correct direction, 0 calls on the null comparison.

```python
import pandas as pd
import numpy as np

se = pd.read_csv('rmats_output/SE.MATS.JC.txt', sep='\t')

def min_per_rep(s):
    return s.str.split(',').apply(lambda x: min(int(v) for v in x))

se['min_inc'] = min_per_rep(se['IJC_SAMPLE_1']).combine(min_per_rep(se['IJC_SAMPLE_2']), min)
se['min_skip'] = min_per_rep(se['SJC_SAMPLE_1']).combine(min_per_rep(se['SJC_SAMPLE_2']), min)

significant = se[
    (se['FDR'] < 0.05) &
    (se['IncLevelDifference'].abs() > 0.10) &
    ((se['min_inc'] + se['min_skip']) >= 10)
].copy()
```

## leafcutter Differential Intron Usage

**Goal:** Detect differential intron-cluster usage annotation-free, capturing novel junctions and complex multi-junction events.

**Approach:** Extract junctions with regtools, cluster introns by shared splice sites, run cluster-level Dirichlet-multinomial test.

```bash
# Strand comes from the XS tag: the BAMs need it (STAR --outSAMstrandField intronMotif); without XS regtools writes strand '?' and leafcutter drops every junction
for bam in *.bam; do
    regtools junctions extract -a 8 -m 50 -s XS "$bam" -o "${bam%.bam}.junc"
done
ls *.junc > juncfiles.txt

# LEAFCUTTER = a clone of github.com/davidaknowles/leafcutter (the R package does not install these scripts on PATH)
python $LEAFCUTTER/clustering/leafcutter_cluster_regtools.py \
    -j juncfiles.txt \
    -o leafcutter \
    -m 50 \
    -l 500000
```

- `-m 50` (min reads per cluster) yields 0 clusters on shallow data (a ~100k-read chrX subset needed `-m 10`).
- Contigs other than `chr1..22,X,Y` / `1..22,X,Y` (e.g. a custom or scaffold reference) are dropped silently and give 0 introns: add `-k True`.
- The counts-table column names are the `.junc` file names without `.junc`.

```bash
# groups.txt: two columns, no header; column 1 must equal the counts-table column names (the .junc basenames, e.g. G1_rep1), column 2 the group
# Extra columns = confounders (see Confounder Handling)
Rscript $LEAFCUTTER/scripts/leafcutter_ds.R --num_threads 4 -i 3 -g 3 -c 10 \
    --exon_file gencode_exons.txt.gz \
    leafcutter_perind_numers.counts.gz groups.txt -o ds_results
```

```r
cluster_sig <- read.table('ds_results_cluster_significance.txt', header = TRUE, sep = '\t')
intron_effects <- read.table('ds_results_effect_sizes.txt', header = TRUE, sep = '\t')

sig_clusters <- subset(cluster_sig, p.adjust < 0.05)
```

**Sample-size flags.** The defaults `-i 5 -g 3 -c 20` stop with `The number of samples in the smallest group is less than min_samples_per_intron` whenever a group has fewer than 5 samples (they stopped at 2v2, 3v3 and 4v4; 5v5 ran). Set `-i` and `-g` to the smaller group size and `-c` (reads per sample) to about 10:

| Design | Flags | Planted 24 strong / 6 weak events found | Null (A vs A') |
|--------|-------|------------------------------------------|----------------|
| 2v2 | `-i 2 -g 2 -c 10` | 21 / 1 | 0 calls |
| 3v3 | `-i 3 -g 3 -c 10` | 22 / 3 | 0 calls |
| 4v4 | `-i 4 -g 4 -c 10` | 23 / 2 | 0 calls |
| 5v5 | `-i 5 -g 5 -c 10` | 24 / 4 | 0 calls |

Simulated data, p.adjust < 0.05 and max |ΔPSI| > 0.10; 4 of 26 calls at 2v2 and 5 of 30 at 3v3 were false, 1 of 26 at 4v4 and 1 of 29 at 5v5, so treat n<4 calls as needing a second tool. The `--min_samples_per_intron 5` "pre-filter" is the default that stops small designs; do not add it.

`--exon_file` (columns chr, start, end, strand, gene_name) only labels clusters with genes; build it with `Rscript $LEAFCUTTER/scripts/gtf_to_exons.R annotation.gtf.gz gencode_exons.txt.gz` (GTF needs `gene_name` attributes) or omit `-e`. `deltapsi` = group 2 - group 1.

**LeafCutter2** (Buen Abad Najar 2025 *bioRxiv*) extends leafcutter with NMD-aware classification of unproductive splicing — useful when AS-NMD coupling is the question.

## MAJIQ V3 Differential Analysis

**Goal:** Detect differential LSVs with full posterior distributions over ΔPSI; ideal for complex multi-junction events and heterogeneous cohorts.

**Status: not run.** MAJIQ V3 and VOILA are academic-licence downloads (majiq.biociphers.org) and are not installed here. What follows is taken from the MAJIQ 3.0.11.dev7 documentation (biociphers.bitbucket.io/majiq-docs); check every flag with `majiq <command> --help` before use.

- Pipeline (V3): `majiq build <gff3> <experiments_tsv> <output_dir>` writes `splicegraph.zarr` and one `<experiment>.sj` per sample (options include `--min-experiments`, `--mindenovo`, `--simplify`, `--strandness`); `majiq psi-coverage <splicegraph> <psi_coverage> <sj ...>` prepares coverage (`--minreads`, `--minbins`); `majiq deltapsi` (two groups of replicates, `-grp1/-grp2`, `-n NAME1 NAME2`, `--splicegraph`) and `majiq heterogen` (two groups of independent experiments) quantify.
- V3 replaced V2's `.majiq` files, `-c settings.ini` and SQLite splicegraph; V2-era flags (`--minpos`, `--mem-profile`, `-j`) and `voila view` on V3 output are unverified. The documentation says VOILA "currently only supports MAJIQ v2".
- Use HET for n>=10 vs n>=10 cohort designs (clinical, GTEx-style; conservative at n=5-10, where deltapsi reports more) and deltapsi for tightly controlled n=3-5.

## SUPPA2 Differential Analysis

**Goal:** Quick differential splicing from existing transcript quantifications, useful as a sanity check or pilot.

**Approach:** Generate per-condition PSI files from Salmon TPM, then run `diffSplice` with empirical or classical p-values.

```bash
suppa.py generateEvents -i annotation.gtf -o events -f ioe -e SE SS MX RI

for ev in SE A5 A3 MX RI; do
    suppa.py psiPerEvent -i events_${ev}_strict.ioe -e ctrl_tpm.tsv -o ctrl_${ev}
    suppa.py psiPerEvent -i events_${ev}_strict.ioe -e trt_tpm.tsv -o trt_${ev}

    suppa.py diffSplice \
        -m empirical \
        -gc \
        -i events_${ev}_strict.ioe \
        -p ctrl_${ev}.psi trt_${ev}.psi \
        -e ctrl_tpm.tsv trt_tpm.tsv \
        -o diff_${ev}
done
```

The TPM files need a header of sample names only (no `Name`/`transcript` label); generateEvents on an SE-only annotation writes empty A5/A3/MX/RI files. SUPPA2 2.4 imports `multipletests` from a statsmodels module that statsmodels 0.15 removed: use statsmodels <=0.14 (0.14.6 worked).

**Minimum replicates: n>=4 per group, in either mode.** Measured on a simulation (24 strong + 6 weak planted ΔPSI events among 120 genes, TPM from the same molecules; calls at p<0.05 and |ΔPSI|>0.10):

| n per group | empirical: strong found, false calls / all calls, calls on a null comparison | classical: strong found, smallest p |
|-------------|------|------|
| 2 | 19/24, 6/26, 8 | 0/24, p 0.22 |
| 3 | 21/24, 4/28, 1 | 0/24, p 0.077 |
| 4 | 22/24, 1/24, 0 | 21/24, p 0.027 |
| 5 | 21/24, 0/22, 0 | 24/24, p 0.008 |

- **Classical (Wilcoxon) is not a fallback for n<=3.** An unpaired rank-sum test cannot reach p<0.05 there (BH cannot lower it), so it found nothing at 2v2 and 3v3 (real chrX 2v2: smallest p 0.22).
- **Empirical at n<=3** finds most strong events but 14-23% of its calls are false and it calls events on a null comparison (real chrX 2v2, permuted labels: 24 events; true labels: 36). Use it only to screen candidates and confirm with rMATS/leafcutter/Shiba; for n<=3 prefer those tools.
- ΔPSI = group 2 - group 1.

## Shiba for Low-Coverage / Few-Replicate Designs

**Goal:** Detect differential splicing with explicit junction-imbalance correction — addresses a known false-positive source for rMATS-style methods.

**Approach:** Shiba is a pipeline configured via YAML. Install via bioconda, write an experiment table and a config file, then run `shiba.py` (`--mame` = splicing analysis only).

```bash
conda install -c bioconda shiba

# exp.tsv (tab-separated): sample  bam  group  technology (short|long); groups named in config.yaml
# config.yaml (keys used in the checked run):
#   workdir: shiba_out          gtf: annotation.gtf        experiment_table: exp.tsv
#   unannotated: False          minimum_anchor_length: 6   minimum_intron_length: 50
#   maximum_intron_length: 500000   strand: XS            only_psi: False   only_psi_group: False
#   fdr: 0.05                   delta_psi: 0.1             reference_group: Ref   alternative_group: Alt
#   minimum_reads: 10           individual_psi: True       ttest: False      excel: False
shiba.py -p 8 --mame config.yaml
# results: shiba_out/results/splicing/PSI_<SE|FIVE|THREE|MXE|RI|AFE|ALE|MSE>.txt; columns dPSI, q, "Diff events" (Yes/No)
```

- `strand: XS` needs XS-tagged BAMs; regtools 1.0 rejects the numeric `strand: 0/1/2` values Shiba passes through.
- The GTF must contain every event type with reads (Shiba 0.8.2 raised a KeyError on a GTF with SE events only).
- dPSI = alternative - reference group. `snakeshiba.smk` lives in the package's share directory and its config needs a container key, so `shiba.py` is the entry point that ran.
- Shiba (Kubota 2025 *NAR*) is reported to outperform rMATS at n=2 vs n=2 by correcting differential mappability between inclusion and skipping junctions. Not reproduced here: on a simulation without planted junction imbalance it found 22/24 strong events at 4v4 with 2 false calls, no better than rMATS or leafcutter, and made 2 calls on a null 2v2 comparison (separate simulation). Community calibration is still emerging. Full config schema: https://sika-zheng-lab.github.io/Shiba/

## Per-Tool Failure Modes

### leafcutter: Cluster Mis-Topology

**Trigger:** A cluster spans a complex topology (cassette + alternative donor in same cluster).

**Mechanism:** Cluster-level p-value reports "something in this cluster differs" but doesn't indicate which intron drove the change; downstream analysis needs per-intron effect sizes.

**Symptom:** Significant cluster, multiple introns with different ΔPSI directions, ambiguous biological interpretation.

**Fix:** Inspect cluster in leafviz; report per-intron effect sizes from `ds_results_effect_sizes.txt`; map cluster topology to canonical SE/A5SS/A3SS via flanking exon coordinates manually.

## Reconciliation: When Tools Disagree

The two most common short-read tools answer slightly different questions: rMATS classifies on annotated event templates; leafcutter classifies on observed cluster usage. Disagreement is informative.

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| rMATS sig, leafcutter not sig | rMATS junction imbalance OR rMATS event hits annotation that leafcutter clustered differently | Inspect locus in IGV; check Shiba on the same locus |
| leafcutter sig, rMATS not sig | Novel junction not in rMATS annotation; rMATS `--novelSS` may have missed it | Verify `--novelSS` was on; rerun if not |
| Both sig, opposite ΔPSI direction | First check sign conventions: rMATS `IncLevelDifference` = group 1 - group 2; SUPPA2 dPSI and leafcutter `deltapsi` = group 2 - group 1; Shiba dPSI = alternative - reference group (verified on planted data). Then event class mismatch (rMATS calls SE positive, leafcutter sees A5SS shift in same cluster) | Flip one sign; if still opposite, manually map cluster topology to event class |
| Both sig, same direction | High-confidence call | Report; cross-validate with sashimi-plot |
| All tools null but biology suggests change | Underpowered design or wrong regime | Increase replicates; check whether outlier-splicing-detection regime applies |

**Operational rule:** for high-confidence reporting, require concordant detection in two tools from different algorithmic families (event-based + cluster-based, or LSV + isoform-based). Document both calls and any explainable disagreements.

## rMATS Output Columns Reference

| Column | Meaning |
|--------|---------|
| IJC_SAMPLE_1 / SJC_SAMPLE_1 | Comma-delimited inclusion / skipping junction counts per replicate, group 1 |
| IJC_SAMPLE_2 / SJC_SAMPLE_2 | Same for group 2 |
| IncFormLen / SkipFormLen | Effective lengths normalizing PSI for differential mapping opportunity |
| upstreamES/EE, downstreamES/EE | Flanking exon coordinates (genomic order; strand-agnostic in column meaning) |
| exonStart_0base / exonEnd | Cassette exon coordinates (0-based half-open) |
| PValue | LRT p-value of \|ΔPSI\| > cutoff |
| FDR | BH-adjusted PValue within event class |
| IncLevel1, IncLevel2 | Comma-delimited per-replicate PSI values |
| IncLevelDifference | mean(IncLevel1) - mean(IncLevel2); sign matches --b1 - --b2 order |

## Replicate Count and Power

| Design | Expected power for ΔPSI=0.2 |
|--------|------------------------------|
| n=2 vs n=2 | Marginal; many real effects missed; SUPPA2 unusable |
| n=3 vs n=3 | Adequate at moderate coverage; standard |
| n=5 vs n=5 | Good; recommended for publication |
| n=10+ vs n=10+ heterogeneous | MAJIQ-HET designed for this scale |
| Single patient vs n=20+ controls | Outlier regime (leafcutterMD, FRASER2); see outlier-splicing-detection |

For an effect-size of |ΔPSI|=0.10 (typical biological signal), power generally requires n>=4 and >=20 junction reads per replicate. Below this, expect to miss most real changes.

## Significance and Effect-Size Thresholds

| Stringency | \|ΔPSI\| | FDR | Use case |
|------------|----------|-----|----------|
| Lenient | > 0.05 | < 0.10 | Discovery, exploratory, hypothesis generation |
| Standard | > 0.10 | < 0.05 | Publication; default reporting threshold |
| Stringent | > 0.20 | < 0.01 | Validation cohort, follow-up targets |

For MAJIQ: posterior probability `P(|ΔPSI| > 0.2) >= 0.95` is roughly equivalent to standard stringency (a different scale from FDR; confirm by simulation when reporting). Always document tool, threshold, and rationale.

**Biologically meaningful ΔPSI varies by context:**
- A poison exon shift of |ΔPSI|=0.10 can halve functional protein (huge biology, modest number).
- A stoichiometric isoform shift of |ΔPSI|=0.10 may be physiologically silent.
- Therapeutic ASO target: SMA nusinersen aims for ΔPSI~+0.30 in SMN2 exon 7.

## Confounder Handling

**Trigger:** sequencing batch, RIN, library prep date or sex correlates with the comparison. rMATS' default LRT takes no covariates (`--paired-stats` handles pairing only), so hits can be driven by batch: PCA on the PSI matrix shows samples clustering by batch rather than group. **Check this first**; if PC1 separates by batch rather than group, the comparison is confounded.

**Batch identical to (or perfectly nested in) condition cannot be adjusted.** leafcutter_ds.R with `batch = group` still exits 0 and returns significant clusters, with no warning; test the design matrix rank before trusting any covariate model (the snippet below does).

Workarounds for rMATS:
1. **Stratification**: run rMATS within each batch separately and meta-analyze.
2. **Per-event regression with a group term (logit-transformed PSI)**: PSI is bounded [0,1], so logit-transform, fit `logit_psi ~ group + batch (+ RIN)` and test the **group coefficient**. Do not regress out batch alone and rank-test the residuals: with 3 vs 3 the rank-sum test cannot go below p=0.1 (0/20 strong events reached p<0.05), and residualizing on an imbalanced batch removes group signal.
3. **Switch to leafcutter** (R function accepts `confounders` matrix; CLI accepts confounders as additional columns in the groups file).

```python
import numpy as np
import pandas as pd
import patsy
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests

# meta: one row per sample (columns group 0/1, batch, RIN); psi: one row per sample for ONE event (columns psi + the same covariates)
design = patsy.dmatrix('group + C(batch) + RIN', meta)
if np.linalg.matrix_rank(np.asarray(design)) < design.shape[1]:
    raise ValueError('a covariate is aliased with group: this comparison cannot be adjusted')

eps = 1e-3  # logit-transform PSI before modelling (PSI is bounded [0,1])
p = psi['psi'].clip(eps, 1 - eps)
psi['logit_psi'] = np.log(p / (1 - p))
fit = smf.ols('logit_psi ~ group + C(batch) + RIN', data=psi).fit()
p_group = fit.pvalues['group']   # one p per event; then multipletests(p_values, method='fdr_bh') across events
```

Power is limited by residual degrees of freedom (samples minus model terms): on a simulation with an imbalanced batch plus RIN, 3 vs 3 found 0/21 strong events at BH q<0.05 (9/21 at nominal p<0.05; residual df=2), 5 vs 5 found 21/21 with 1 false call among 46 null events. With few replicates drop RIN, or use stratification.

**leafcutter** accepts confounders two ways (verified in the script source and by a run):
- **R function**: `differential_splicing(counts, x, confounders=numeric_matrix)` accepts a numeric covariate matrix
- **CLI script**: `leafcutter_ds.R` reads confounders from **additional columns in the groups file** (3rd, 4th, ... columns), NOT from a `--confounders` flag; use labels like `b1`/`b2` for categorical variables (numeric columns are treated as continuous)

**MAJIQ** does not accept arbitrary confounders (its documentation lists a separate `moccasin` batch-correction step; unverified); use stratification or switch tool.

## Multi-Group / Multi-Factor Designs

| Design | Approach |
|--------|----------|
| 3 groups (e.g. drug A, drug B, control) | Pairwise rMATS or leafcutter; OR limma/DESeq2 on logit-PSI matrix |
| Time-course (e.g. 0h, 6h, 24h) | DEXSeq on event counts with time as factor; or limma::lmFit on PSI matrix |
| 2x2 factorial (genotype × treatment) | DEXSeq with interaction term; rMATS pairwise on interaction subsets |
| Continuous covariate (dose, age) | limma::lmFit on logit-PSI ~ covariate |

For complex designs, custom regression on the PSI matrix is more flexible than rMATS/leafcutter pairwise.

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| `rMATS: numpy.AxisError` | rMATS version mismatch with numpy >=2.0 | Pin numpy<2.0 or update rMATS-turbo to >=4.3 |
| rMATS exits 0, `SE.MATS.JC.txt` has only a header (or no FDR column) | `-t` does not match the reads, `--readLength` mismatch without `--variable-read-length`, or `--paired-stats` without PAIRADISE (`no package called PAIRADISE` in the log) | Fix the flag / install PAIRADISE; assert rows + FDR after every run |
| `leafcutter_ds.R: The number of samples in the smallest group is less than min_samples_per_intron` | Default `-i 5` with a group of <5 samples | `-i N -g N` with N = smaller group size (leafcutter section) |
| `leafcutter_ds.R: undefined columns selected` | Sample names in groups.txt do not equal the counts-table column names (.junc basenames) | Rename groups.txt entries to the .junc basenames |
| leafcutter clustering writes 0 introns | Contig not `chr1..22,X,Y`/`1..22,X,Y`; `-m` too high for the depth; or `.junc` strand `?` (no XS tag) | `-k True`; lower `-m`; use XS-tagged BAMs |
| SUPPA2 `ImportError` on `multipletests` | statsmodels 0.15 | Install statsmodels <=0.14 |
| `SUPPA2: no events with sufficient coverage` | Salmon/kallisto TPM filter too strict upstream | Lower upstream TPM threshold; verify event annotations |
| `regtools: too many open files` | Many BAMs in one batch | `ulimit -n 4096` or batch in groups |

## Result Prioritization

**Goal:** Rank events by combined statistical and biological significance for follow-up.

**Approach:** Composite score combining FDR and effect size, then enrich for biology (RBP binding, NMD sensitivity, conservation, disease relevance).

```python
import pandas as pd
import numpy as np

# `significant` = the filtered rMATS table from the rMATS section
sig = significant.copy()
sig['score'] = -np.log10(sig['FDR']) * sig['IncLevelDifference'].abs()
sig['exon_length'] = sig['exonEnd'] - sig['exonStart_0base']
sig['nmd_likely'] = (sig['exon_length'] % 3 != 0)
top_events = sig.nlargest(50, 'score')
```

Cross-reference top hits with:
- **eCLIP/ENCODE RBP target databases** (POSTAR3, oRNAment, RBP2GO) -> candidate trans-regulators
- **Disease-specific signatures**: SF3B1 cancer signal = cryptic 3'ss ~10-30nt upstream of the canonical one (MDS/CLL/UM); TDP-43 cryptic exons (UNC13A, STMN2, ATG4B) for ALS/FTD
- **Conservation**: VastDB cross-species PSI for evolutionary support
- **Splice-site predictions**: SpliceAI scores for the involved sites (see splice-variant-prediction)

## Common Pitfalls

- **Junction read imbalance** (cassette exon flanks have unequal mapping opportunity) inflates rMATS false positives; Shiba explicitly corrects this.
- **Forgetting NMD direction** — increased PSI of a poison exon decreases protein. Always check whether the alternative form is PTC-introducing using ORF-aware annotation.
- **Cryptic splicing in TDP-43 loss / SF3B1-mutant samples** — annotation-bound tools (rMATS, SUPPA2) miss these; need leafcutter or MAJIQ with denovo mode.
- **Forgetting strand** — wrong `--libType` halves usable junctions. Confirm with RSeQC `infer_experiment.py`.

## Related Skills

- splicing-quantification - PSI estimation per event; foundational
- splicing-qc - Run BEFORE differential to verify library, depth, strandedness; avoid downstream surprises
- isoform-switching - DTU framework with NMD/ORF/domain consequences; complementary to event-level
- sashimi-plots - Visualize differential events for QC and reporting
- outlier-splicing-detection - Single-sample-vs-cohort regime (FRASER2/DROP); use when not 2-group
- splice-variant-prediction - SpliceAI / Pangolin for variant-driven mechanistic explanation of differential events
- long-read-splicing - Differential analysis from full-length isoforms; use when short-read insufficient
- read-alignment/star-alignment - STAR 2-pass cohort-style required upstream

## References

- Shen et al 2014 *PNAS* - rMATS original
- Wang et al 2024 *Nat Protoc* - rMATS-turbo
- Li et al 2018 *Nat Genet* - leafcutter (Dirichlet-multinomial GLM)
- Buen Abad Najar et al 2025 *bioRxiv* - LeafCutter2 (NMD-aware unproductive splicing)
- Vaquero-Garcia et al 2016 *eLife* - MAJIQ LSV framework
- Vaquero-Garcia et al 2023 *Nat Commun* - MAJIQ-HET heterogeneity module
- Aicher, Slaff, Jewell, Barash 2024 *bioRxiv* - MAJIQ V3
- Trincado et al 2018 *Genome Biol* - SUPPA2
- Kubota et al 2025 *NAR* - Shiba (junction-imbalance correction)
- Olofsson et al 2023 *Biochem Biophys Res Commun* 653:31-37 - benchmark across tools
- Tran et al 2025 *WIREs RNA* - methodology review
- Brown et al 2022 *Nature* - UNC13A cryptic exon (TDP-43 / ALS)
- Klim et al 2019 *Nat Neurosci* - STMN2 cryptic splicing (ALS)
- Darman et al 2015 *Cell Rep* - SF3B1 cryptic 3'ss
