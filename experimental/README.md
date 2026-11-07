# Experimental normalized Skills

These 24 `ca-` (context-aware) Skills preserve partially complete revisions of
the audited collection. Their short entry points route to task instructions,
scripts, and deeper references. They are available for inspection and further
work; publication does not certify their scientific behavior.

One has a historical candidate-ready relay result, three have limited evaluation
results without an audit, and twenty were normalized without an audit. Read
each linked status note for executed surfaces, known defects, and missing work.
Do not transfer the original Skill's audit score to a changed candidate.

Sam authorized publication of these drafts and resumed normalization-only work
on the other 171 optimized shelf Skills on 2026-11-07. Scientific audits and
Open Science Marketplace follow-up remain outside this resumed scope.

The canonical drafts live in optimized-scientific-skills; bioSkills-Improved
receives matching copies. Existing published versions remain available under
their original IDs. Install a candidate explicitly by copying its entire
`skills/ca-<id>` folder into your agent's skills directory. Each folder's
frontmatter name matches its `ca-` directory name.

| Candidate | Recorded status |
| --- | --- |
| [ca-bio-workflows-crispr-screen-pipeline](skills/ca-bio-workflows-crispr-screen-pipeline/STATUS.md) | Through the full relay; candidate-ready in `audits/skills/bio-workflows-crispr-screen-pipeline/candidate@e609525649bd-reaudit-delta-001` (final 89). Only the frontmatter `name` differs from the audited bytes. |
| [ca-bio-differential-expression-deseq2-basics](skills/ca-bio-differential-expression-deseq2-basics/STATUS.md) | Not audited. Eval: pseudobulk task 5 of 5 with a cheap model (old layout 2 of 5); routing check 7 of 7 routes. Eleven older reference files only lightly cleaned. |
| [ca-bio-proteomics-differential-abundance](skills/ca-bio-proteomics-differential-abundance/STATUS.md) | Not audited. Eval: 5 of 5 (old layout 2 of 5). DEqMS, proDA, msqrob2 and MSstats scripts and `shrink=ashr` never executed. No routing cases. |
| [ca-bio-machine-learning-model-validation](skills/ca-bio-machine-learning-model-validation/STATUS.md) | Not audited. Eval: 3 of 5 (old layout 2 of 5); routing check 9 of 9 routes. Known defects: the decision-curve route needs a predictions file nothing produces; one route overstates that leave-one-out cannot give an AUC. |
| [ca-bio-duplicate-handling](skills/ca-bio-duplicate-handling/STATUS.md) | Normalized only, not audited. Executed: all five scripts (samtools, pysam, umi_tools, fgbio on staged BAMs). Unexecuted: `umi_tools_dedup.sh scrna` success path, alternative-marker commands. |
| [ca-bio-variant-calling-filtering-best-practices](skills/ca-bio-variant-calling-filtering-best-practices/STATUS.md) | Normalized only, not audited. All scripts unexecuted (`hard_filter.sh`, `genotype_filter.sh`, `somatic_filter.sh`, `cyvcf2_filter.py`, `filter_stats.sh`); VQSR commands shown only. |
| [ca-bio-single-cell-clustering](skills/ca-bio-single-cell-clustering/STATUS.md) | Normalized only, not audited. All scripts unexecuted (`cluster_scanpy.py`, `cluster_seurat.R`, `validate_partition.py`, `subcluster.py`). |
| [ca-bio-single-cell-preprocessing](skills/ca-bio-single-cell-preprocessing/STATUS.md) | Normalized only, not audited. All scripts unexecuted (`preprocess_scanpy.py`, `preprocess_seurat.R`, `soupx.R`). |
| [ca-bio-single-cell-data-io](skills/ca-bio-single-cell-data-io/STATUS.md) | Normalized only, not audited. All scripts unexecuted (`load_10x.py`, `load_10x.R`, `h5ad_inventory.py`). |
| [ca-bio-proteomics-proteomics-qc](skills/ca-bio-proteomics-proteomics-qc/STATUS.md) | Normalized only, not audited. Executed once earlier: `raw_qc.py`, `matrix_metrics.py`, `pca_batch.py`, `tmt_balance.py`, `diann_level1.py` (synthetic MaxQuant/TMT, real DIA-NN report). Unexecuted: `ptxqc_report.R`. |
| [ca-bio-analytical-validation](skills/ca-bio-analytical-validation/STATUS.md) | Normalized only, not audited. Executed: `ge_and_poisson.py`, `lod95_probit.py`, `panel_integrated_lod.py` (contrived dilution series; ETAM-seq-derived series). |
| [ca-bio-machine-learning-omics-classifiers](skills/ca-bio-machine-learning-omics-classifiers/STATUS.md) | Normalized only, not audited. Executed once: `logistic_regression.py --signature` on Golub (dense path crashed, fixed, not rerun). Unexecuted: `rf_xgboost_classifier.py`, `batch_checks.py`, `calibration_check.py`; SMOTE path untestable (no imblearn). |
| [ca-bio-variant-annotation](skills/ca-bio-variant-annotation/STATUS.md) | Normalized only, not audited. All scripts unexecuted (`annotate_vep.sh`, `annotate_vcf.sh`, `parse_csq.py`). |
| [ca-bio-machine-learning-prediction-explanation](skills/ca-bio-machine-learning-prediction-explanation/STATUS.md) | Normalized only, not audited. All scripts unexecuted (`shap_tree.py`, `lime_local.py`); lime in no WSL environment. |
| [ca-bio-pathway-go-enrichment](skills/ca-bio-pathway-go-enrichment/STATUS.md) | Normalized only, not audited. All scripts unexecuted (`ora_go.R`, `goseq_ora.R`, `groupgo.R`). |
| [ca-bio-atac-seq-consensus-peakset](skills/ca-bio-atac-seq-consensus-peakset/STATUS.md) | Normalized only, not audited. Scripts: `iterative_overlap.sh` (shelf copy, unchanged), `iterative_overlap.R`, `per_condition_consensus.sh` (new, unexecuted). |
| [ca-bio-atac-seq-motif-deviation](skills/ca-bio-atac-seq-motif-deviation/STATUS.md) | Normalized only, not audited. All scripts unexecuted after the recut (`chromvar_bulk.R`, `chromvar_signac.R`, `chromvar_archr.R`). |
| [ca-bio-atac-seq-atac-qc](skills/ca-bio-atac-seq-atac-qc/STATUS.md) | Normalized only, not audited. No new scripts; shelf copies unchanged (`library_complexity.py`, `atac_qc_metrics.R`, `encode_tss_enrichment.py`, `aggregate_qc.py`). |
| [ca-bio-single-cell-multimodal-integration](skills/ca-bio-single-cell-multimodal-integration/STATUS.md) | Normalized only, not audited. New, unexecuted: `cite_seq_wnn.R`, `cite_seq_muon.py`, `multiome_wnn.R`, `mofa_muon.py`, `glue_unpaired.py`; shelf copies unchanged: `totalvi_cite_seq.py`, `multivi_mosaic.py`, `seurat_bridge_integration.R`. |
| [ca-bio-alignment-sorting](skills/ca-bio-alignment-sorting/STATUS.md) | Normalized only, not audited. New, unexecuted: `sort_bam.sh`, `merge_bams.sh`, `check_order.py`; `align_sort.sh` is the shelf `sort_pipeline.sh`, renamed. |
| [ca-bio-vcf-statistics](skills/ca-bio-vcf-statistics/STATUS.md) | Normalized only, not audited. All scripts unexecuted (`vcf_qc.sh`, `vcf_stats.py`, `missing_hwe.sh`, `het_allele_balance.sh`, `identity_qc.sh`, `stratified_stats.sh`). |
| [ca-bio-variant-normalization](skills/ca-bio-variant-normalization/STATUS.md) | Normalized only, not audited. All scripts unexecuted (`normalize_vcf.sh`, `check_normalization.py`, `check_ref.sh`, `compare_callers.sh`). |
| [ca-bio-microbiome-amplicon-processing](skills/ca-bio-microbiome-amplicon-processing/STATUS.md) | Normalized only, not audited. All scripts unexecuted (`remove_primers.sh`, `dada2_asv.R`, `decontam_filter.R`); `dada2_asv.R --single` is new code. |
| [ca-bio-microbiome-taxonomy-assignment](skills/ca-bio-microbiome-taxonomy-assignment/STATUS.md) | Normalized only, not audited. All scripts unexecuted (`dada2_assign.R`, `idtaxa_classify.R`, `qiime2_region_nb.sh`, `qiime2_vsearch.sh`, `filter_organelle.sh`). |

The [manifest](CANDIDATES.json) records source commits, per-file hashes,
upstream paths and status. GPTomics/bioSkills supplied the original scientific
Skills under the repository [MIT license](../LICENSE). The two historical
DESeq2 pointer edits remain unapplied to the published collection; they are
preserved as [shelf-pointer-edits.patch](https://github.com/mrsonord2240/optimizing-agent-science-skills/blob/e9dcbbb67b1e3f02922e5d84b6778d757dce08d1/skills/shelf-pointer-edits.patch).
