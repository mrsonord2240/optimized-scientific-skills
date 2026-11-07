# Candidate synthetic tasks for the Skills eval (draft, 2026-09-19)

Goal: find out whether the optimized Skills have utility over upstream bioSkills and over nothing. Not replicating benchmarks.bio.

Nothing below is built. Pick, cut, or reshape before I build anything.

## Rules every task follows

- **Synthetic, offline, ground truth by construction.** We generate the data from a known model, so the answer key is exact and nothing needs the internet.
- **Trap from the domain, not from our Skill text.** Each task names its published pitfall. That keeps the test independent of how we edited the Skills.
- **Deterministic grader.** Numeric tolerance, exact label set, or precision/recall on a gene/feature list (LatchBio's graders, unchanged). Where none fits I add a small custom check and say so.
- **Reference solvers included.** One correct, plus one solver per trap that must fail. This proves the key and grader work; it does not prove an LLM finds the task hard.
- **Headroom check before it counts.** A task is kept only if the no-skill baseline for the chosen model is neither ~0% nor ~100%. Tasks at either end are dropped from the headline.

Two kinds of task, always reported separately:
- **Independent (I-)**: built from a literature pitfall.
- **Regression (R-)**: built so that following the *upstream* Skill's known defect gives a wrong answer, taken from our fix logs. These favour our version by construction, so they measure "did our edits matter", not general utility.

## Independent tasks

Effort: S = a day-equivalent of my work or less, M = a few times that, L = large.

| id | Skill area (ids) | scenario and data to simulate | pitfall (source) | answer / grader | effort |
| --- | --- | --- | --- | --- | --- |
| I1 | machine-learning-model-validation | 60 samples x 5,000 noise features (true AUC 0.5) plus a second set with 5 real features and repeated samples per patient. Report the honest generalisation AUC of "select top-50 features, then classifier". | feature selection outside the CV loop; patient-level leakage (Ambroise & McLachlan 2002) | two numeric AUCs, each within a band; numeric_tolerance | S |
| I2 | differential-expression-deseq2-basics, single-cell-* | Multi-donor single-cell counts, imbalanced cells per donor and condition, known DEG set. Report DEGs. | cells treated as replicates (Squair 2021) | gene list vs truth; marker_gene_precision_recall, or a count band | S |
| I3 | pathway-go-enrichment, pathway-gsea, pathway-* | DE list from a tissue where only ~4,000 genes are detected, a custom GMT with 50 pathways, one pathway enriched only because of the wrong background. | ORA against all genes instead of the expressed universe (Timmons 2015) | significant-pathway set, exact; label_set_jaccard | S |
| I4 | experimental-design-power-analysis | Negative-binomial parameters plus a genome-wide (or proteome-wide) multiplicity rule stated in the task. Ask for n per group at 80% power. Truth by simulation. | multiplicity ignored, understating n several-fold; t-test power on counts | numeric n, +/-20%; convention risk, so the task states every convention | M |
| I5 | crispr-screens-copy-number-correction, mageck, hit-calling | 500 genes x 4 guides, 40 essential genes, one high-copy amplicon of 30 non-essential genes that deplete anyway. Copy-number segments provided. Single cell line. | copy-number-driven false dependencies (Aguirre 2016, Meyers 2017) | hit list; precision/recall with both thresholds | M |
| I6 | proteomics-differential-abundance, quantification | 1,500 proteins, 6 vs 6, MNAR missingness at low abundance, per-sample loading differences, 10% truly changed. | min/zero imputation false positives; sample loading not centred | protein list; precision/recall | M |
| I7 | metabolomics-isotope-tracing | Simulated mass-isotopologue distributions for a C6 metabolite with known true labelling, natural abundance and tracer impurity convolved in. Report the corrected fractions. | uncorrected natural abundance (Buescher 2015) | corrected M+n fractions, +/-0.02; numeric_tolerance | S |
| I8 | microbiome-differential-abundance | Taxa counts as a phyloseq `.rds`, cases with 3x total microbial load, 5 truly differential taxa, some zero-variance taxa. | relative-abundance tests flag compositional artefacts (Nearing 2022) | taxa list; precision/recall | M |
| I9 | variant-normalization, vcf-manipulation | Two callers' VCFs describing the same variants differently: unnormalised indels, multi-allelic sites, split MNVs. Report the concordant and caller-specific counts. | comparing raw positions and alleles | three counts, exact; numeric_tolerance | S |
| I10 | causal-genomics-colocalization-analysis | Reuse the MR region simulator: four regions, two traits each: shared causal SNP; two distinct SNPs in high LD; two independent signals sharing one; no overlap. Call colocalisation per region. | single-causal-variant assumption gives false colocalisation | four labels, exact; dict match | M |
| I11 | molecular-standardization, similarity/substructure | 300 SMILES derived from 100 parents by adding counter-ions, changing charge state, reordering atoms. Report unique parent count. Needs rdkit. | canonical SMILES alone does not merge salts or charge forms | count, +/-2; task states the policy | S |
| I12 | qsar-modeling | 600 molecules from 12 analogue series, activity from series plus substituent terms. Report honest RMSE for new series. Needs rdkit. | random-split CV leaks analogues (scaffold or group split needed) | numeric RMSE band | M |
| I13 | population-genetics-rare-variant-association | Gene-level rare-variant data, some genes with opposing-direction effects, stratified ancestry. Which genes are significant? | burden test loses power to opposing effects; stratification | gene set; precision/recall | M |
| I14 | experimental-design-batch-design / randomization-blocking | A confounded 24-sample layout across 3 batches. Return a corrected allocation. | condition fully confounded with batch | custom checker on the allocation (no stock grader fits) | S |

Deliberately left out: phylogenetics (heavy tools, needs dating software), docking/pose validation (3D, heavy), single-cell doublets and integration (scanpy stack, later).

**Cannot be tested offline** (they call databases or download data): entrez-*, geo-data, sra-data, biomart-queries, ncbi-datasets-cli, batch-downloads, clinical-databases-*, gnomad, clinvar, dbsnp. About 12 Skills. I'd exclude them and say so.

## Regression tasks (from our fix logs)

Each is an existing task above with a specific input that trips the upstream defect.

| id | on | upstream defect our fix corrected (fix log) | what the task needs so the defect matters |
| --- | --- | --- | --- |
| R1 | I6 | proteomics-differential-abundance: example "median centering" divided log values instead of subtracting; `treat()` dropped trend and robust | strong per-sample loading differences |
| R2 | I4 | experimental-design-power-analysis: proteomics route dropped proteome-wide multiplicity, understating n 3-7x | proteome-wide test count in the scenario |
| R3 | I9 | variant-normalization: atomize-before-split emitted spurious `*` records; csq needed `--phase` | multi-allelic and overlapping-deletion records |
| R4 | new: labelling time-course | metabolomics-isotope-tracing: steady state judged from only the last two timepoints | a still-decelerating series; ask "is steady state reached at 24 h?" |
| R5 | I5 | crispr-screens-copy-number-correction: decision tree sent single-cell-line screens to Chronos, which needs at least 3 lines | one cell line only |
| R6 | I8 | microbiome-differential-abundance: LinDA crashed on a real phyloseq object; ZicoSeq crashed on zero-variance features | phyloseq `.rds` input, zero-variance taxa |
| R7 | I10 | colocalization: shipped `coloc.susie` examples crashed; PP.H3 trigger overstated | a region needing the SuSiE route |

R1, R3, R5, R6 and R7 are partly *runnability* differences. The upstream recipe throws an error and the agent must recover. That is a real utility difference, but a different one from "gives a wrong answer".

## Sandbox additions

- Bioconductor and CRAN binaries through r2u apt: limma, edgeR, DESeq2, phyloseq, ALDEx2, MicrobiomeStat, SKAT (if available there).
- pip: rdkit, mageck, isocor.
- Already present: MendelianRandomization, susieR, coloc, TwoSampleMR, plink, bcftools.
- Image rebuild is about the same 15-20 minutes as before.

## Suggested order (first batch of 8)

I1, I2, I3, I7, I9 are small and cheap to build. Then I5, I6, I8 for the larger areas. That covers ML, single-cell, pathways, metabolomics, variants, CRISPR, proteomics and microbiome, and includes R1, R3, R5, R6 embedded. I10 (colocalization) reuses the existing simulator and would be next.

Known risks:
- I1, I2 and I3 are the ones a mid-tier model may already do well, so they are the likeliest to sit at ceiling with no skills.
- The Skills' descriptions cost roughly 10-12k tokens per step in the Skill arms; more tasks means more spend, but at deepseek-flash prices still a few dollars per full pass.
