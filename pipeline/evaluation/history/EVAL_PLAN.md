# Skills eval plan: optimized vs upstream vs none (2026-09-19)

Question: do the refined Skills in `optimized-scientific-skills` beat the original `GPTomics/bioSkills`, and do either beat no Skills, on benchmarks.bio problems?

## What benchmarks.bio is

LatchBio's verifiable agentic benchmarks. Each problem = a data snapshot + natural-language task + deterministic grader on a structured `<EVAL_ANSWER>` JSON. Cloned here (public subsets only):

| repo | full size | public here |
| --- | --- | --- |
| `scbench` (scRNA-seq) | 195 | 6 |
| `spatialbench` | 159 | 15 |
| `variantbench` | 118 | 8 |
| `epibench` | 106 | 7 |
| `latch-eval-tools` | graders + `EvalRunner` + agent envs | n/a |

The full sets are withheld to prevent contamination. **The usable pool is the public subset only.**

## Coverage gap

The benchmarks touch single-cell, variant/VCF, and statistical-genetics Skills. They do not touch proteomics, metabolomics, CRISPR screens, microbiome, phylogenetics, pathway, or chemoinformatics, which is most of the 126 Skills. A win here says nothing about those.

## Task selection

Arm B (upstream) is the same skill ids at `GPTomics/bioSkills@d91ed3d`. Where our Skill is `licence-declaration-only`, arms A and B are byte-identical apart from `license: MIT`, so those tasks are free A/A noise controls.

### Tier 1: Skill modified by us, so A can differ from B

| eval | repo | Skills it should trigger | trap it tests |
| --- | --- | --- | --- |
| `fine_mapping_effector_gene_prioritization` | variantbench | fine-mapping, effector-gene-prioritization, variant-annotation | benign-predicted missense in the credible set; no cis-eQTL. Nearest-gene / cis-eQTL / top-trans heuristics all fail |
| `cis_mr_bcaa` | variantbench | mendelian-randomization, proteome-mr-drug-target (unmodified), colocalization | cis-MR for a drug target (BCKDK -> T2D) |
| `small_variant_calling_mnv_codon_consequence` | variantbench | variant-annotation, vcf-manipulation, variant-normalization, variant-calling-filtering | per-variant annotation misses phased same-codon MNV consequences |
| `chry_father_son_concordance_11` | variantbench | vcf-basics, vcf-statistics, vcf-manipulation | **not yet read closely; confirm fit before committing** |
| `sv_chromosome_end_signal` | variantbench | vcf-manipulation, vcf-statistics | **not yet read closely; confirm fit before committing** |

### Tier 2: Skill unmodified, so A/A controls and Skills-vs-none only

| eval | repo | Skills | trap |
| --- | --- | --- | --- |
| `DE01_pseudobulk_de` | scbench | deseq2-basics | single-cell t-test/Wilcoxon instead of donor-aware pseudobulk |
| `NRM01_sparse_normalization` | scbench | single-cell-preprocessing | LogNormalize on shallow data. Its own notes say low discrimination |
| `dr_05_pca_preprocessing_sentinels` | scbench | single-cell-preprocessing | preprocessing sentinels |
| `bd_rhapsody_tnbc_panel_aware_qc` | scbench | preprocessing, doublet-detection | panel-aware QC |
| `T04a_endothelin_niche_sources` | scbench | clustering, markers-annotation | exploratory clustering to cell typing |
| `qc_01_endothelial_doublet_assessment` | spatialbench | doublet-detection | doublets vs real co-expression |
| `DR01_batch_confounded_clustering`, `NORM01_batch_correction`, `curio_ovary_batch_driven_clustering` | spatialbench | batch-integration, clustering | batch-driven structure |

### Excluded

All of `epibench` (bulk ATAC/ChIP/CUT&Tag; we ship no Skill for it), `tapestri_ccus_clustering` (scDNA), `ont_rfc1` (repeat expansions), both neoantigen evals, the xenium niche tasks.

Result: about 15 tasks, 5 of them where A and B can genuinely differ.

## The three arms

| arm | workspace `skills/` | prompt addition |
| --- | --- | --- |
| C0 none | absent | none |
| C1 upstream | 126 ids from `bioSkills@d91ed3d` | identical stub, listing name + description |
| C2 optimized | 126 ids from `optimized-scientific-skills@868f7d3` (pin the commit) | identical stub |

Everything else is held constant: model, temperature, step and token budget, container image, data, grader.

Deliver all 126 in both arms, not only the "right" Skill, so the test includes the agent choosing which Skill to open.

Stub, verbatim in both skill arms: "Reference Skills are in ./skills/<name>/SKILL.md. Read any that are relevant before you start." Followed by the name and description lines.

## Harness

Primary: **mini-swe-agent** (the benchmark's own default, so results sit next to their published table). It goes through litellm, so `openrouter/<model>` with `OPENROUTER_API_KEY` works. Bash-only agent, which is why the Skills are plain files in the workspace.

Secondary, only if progressive-disclosure triggering matters: Claude Code pointed at OpenRouter's Anthropic-compatible endpoint, Skills in `.claude/skills/`. Verify that route before relying on it.

Run inside Docker or WSL. The tasks need bcftools, samtools, scanpy, susieR, and so on, and the same image must serve all three arms. `latch-eval-tools/agent_env/` is a starting point.

## Model

A mid-tier model whose **no-skill baseline lands at roughly 30-60%** on these tasks. Below that there is no headroom to see Skills help; above it, a ceiling. For reference, the published mini-swe-agent scores on scbench are around 50% for a Sonnet-class model and about 24% for the weakest. Run C0 on all tasks first with two candidate models, and pick on the baseline.

## Protocol

1. Strip `notes`, `canary`, and `data_node` from what the agent sees. `notes` contains the answer and the trap.
2. No internet in the container.
3. 5 repeats per task per arm: about 15 x 3 x 5 = 225 runs. Published cost columns are about $0.3-1.5 per task-run for mid-tier models; Skill context adds to that. Pilot with 1 repeat before committing.
4. Grade with `latch-eval-tools` graders, unmodified.
5. Log every trajectory. Compute **uptake**: did the agent `cat` a `skills/*/SKILL.md`, and which one. If uptake is near zero, the arms are the same and no result means anything.

## Analysis (pre-register before running)

- Primary: mean pass rate across tasks, C2 vs C1 vs C0, paired by task. Bootstrap over tasks, with within-task repeat resampling.
- Noise floor: C1 vs C2 on the Tier 2 tasks. They should differ by about 0. The spread there is the smallest effect you can trust on Tier 1.
- With 5 Tier 1 tasks, only large effects (15-20 points or more) are detectable. Report per-task pass rates, not only the mean.
- Secondary: pass rate split by uptake vs no uptake. Selection-biased, so read as descriptive.
- Drop tasks at 0% or 100% in all three arms before computing the headline.

## Risks

- **Data access.** Every eval points at `latch://<id>.node`; `epibench/scripts/download_manifest_data.py` shells out to `latch cp`. Needs the Latch CLI and, very likely, a Latch login. **Not yet verified that the data is publicly downloadable.** This gates everything.
- **Contamination.** The public evals, `notes` and trajectories are on GitHub. Equal across arms, so it shifts the level rather than the delta, but it shrinks headroom.
- **Small n.** About 15 tasks; do not over-read.
- **Task fit.** `chry_*` and `sv_*` not read closely.
