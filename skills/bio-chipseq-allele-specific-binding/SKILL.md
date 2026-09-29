---
name: bio-chipseq-allele-specific-binding
category: Data Analysis
description: Detects allele-specific transcription factor or histone modification binding from heterozygous-variant ChIP-seq using WASP (reference-bias filter), RASQUAL (joint QTL and allelic testing), BaalChIP (Bayesian beta-binomial testing with measured RAF or gDNA correction), and an explicitly external-only AlleleSeq route. Handles imprinted-locus awareness, X-inactivation artifacts, cancer allele-dose imbalance, and downstream caQTL or bQTL integration. Use when identifying variants with allelic effects on TF binding, fine-mapping causal regulatory variants, validating sequence-model predictions, or characterizing cis-acting regulatory effects.
tool_type: mixed
primary_tool: WASP
license: MIT
author: GPTomics
---

## Version compatibility

The executable contracts were checked against WASP v0.3.4 (`f980683`), RASQUAL commit `5aa553c`, and BaalChIP 1.36.0 from Bioconductor 3.22 on R 4.5. The AlleleSeq2 route is documentation-only against [`trgaleev/AlleleSeq2`](https://github.com/trgaleev/AlleleSeq2) commit `cfe8acf` until its official legacy prerequisites are supplied. The prepared environment used samtools/bcftools 1.21 and pysam 0.22.1.

Before running a workflow, inspect the installed version and its help or API documentation. In particular, verify BaalChIP argument names against `?getASB`; releases may differ on `Iter` versus `nIter`.

# Allele-specific binding from ChIP-seq

Compare reference- and alternate-allele ChIP-seq read counts at heterozygous variants to identify cis-acting differences in transcription-factor or histone-mark binding.

Never interpret raw allelic imbalance without addressing these confounders:

1. Reference-allele mapping bias. Apply WASP before count-based testing, or use a method such as RASQUAL that models mapping bias explicitly.
2. Imprinted loci. Flag or filter loci such as H19, IGF2, MEG3, MEG8, and KCNQ1OT1 because constitutive allele skew is not differential binding.
3. X-inactivation. Exclude chrX in female samples from the ordinary autosomal analysis, or analyze it separately with an appropriate model.
4. Copy-number imbalance. In cancer samples, use a copy-number-aware analysis or exclude copy-number-altered regions.

## Route to the needed resource

- Read [methods-and-commands.md](references/methods-and-commands.md) for installation clues, the complete WASP, RASQUAL, and AlleleSeq command surfaces, input formats, and method-selection detail.
- Read [pitfalls-and-troubleshooting.md](references/pitfalls-and-troubleshooting.md) when diagnosing reference skew, chromosome mismatches, convergence failures, sparse coverage, X-inactivation, imprinting, or copy-number artifacts.
- Run or adapt [baalchip_workflow.R](scripts/baalchip_workflow.R) for the complete cancer-oriented BaalChIP example. It retains the provider's executable example as a script instead of duplicating it inline.

## Required inputs and preflight

Gather and record:

- coordinate-sorted ChIP-seq BAMs and their indexes;
- peak calls for the same reference assembly;
- a sample-specific heterozygous-variant VCF or table, preferably from matched-normal calling for cancer samples;
- sample sex and the intended treatment of chrX;
- an imprinted-locus interval set mapped to the same assembly;
- a genome-build-matched blacklist;
- either measured per-variant reference allele frequencies (`RAF`) or matching gDNA BAMs for BaalChIP cancer correction; ordinary population allele frequency (`AF`) and total/minor-copy-number segments are not substitutes;
- phased genotypes when using RASQUAL population analyses or AlleleSeq personalized genomes.

Confirm that chromosome names, coordinates, declared assembly, and sample identifiers agree across BAM headers, peaks, variants, blacklist, imprinted loci, and the selected RAF/gDNA source. Do not silently substitute a population SNP panel for the sample's own heterozygous variants. Fail closed when correction inputs are missing: setting `RAFcorrection=TRUE` without real RAF or gDNA is not evidence of copy-number correction.

## Choose the analysis path

| Situation | Path | Reason |
|---|---|---|
| Single sample or replicate group with measured RAF or matching gDNA | WASP-filtered BAM -> BaalChIP | Bayesian beta-binomial inference with an explicit supported allele-dose correction input |
| Population cohort with phased genotypes and cis-QTL goals | RASQUAL | Joint association and allelic signal with per-feature bias parameter `phi` |
| One sample with high-quality phased genotype and maximum mapping-bias control | AlleleSeq | Aligns to personalized maternal and paternal genomes |
| Simple count-table check after valid bias correction | Beta-binomial or chi-squared test | Useful as a transparent secondary analysis, not a replacement for bias correction |

WASP is preprocessing, not the ASB test. It swaps alleles in reads overlapping heterozygous SNPs, remaps them, and discards reads that do not return consistently. Mapping-bias filtering can materially reduce usable depth; measure and report actual stage counts rather than assuming a fixed attrition percentage. RASQUAL's `phi` is the documented alternative that models bias within the joint test.

## Core workflow

### 1. Establish heterozygous variants

Use a sample-specific callset, restrict to heterozygous small variants, and normalize the representation before deriving tool-specific inputs. For tumor ChIP-seq, prefer variants called from a matched normal rather than inferring germline heterozygosity from copy-number-distorted tumor reads.

### 2. Correct mapping bias

For the default path, run WASP with SNP tables built from the sample-specific VCF:

1. Align reads and coordinate-sort the initial BAM.
2. Identify reads overlapping heterozygous SNPs and generate allele-swapped reads.
3. Remap the swapped reads with the same reference and aligner settings.
4. Keep only reads that map back consistently.
5. Sort and index the final `.wasp.bam` used downstream.

Use the exact command sequence and dependency clues in [methods-and-commands.md](references/methods-and-commands.md). Record input and retained read counts rather than assuming a fixed loss rate.

### 3. Remove or stratify known biological confounders

- Remove imprinted loci before ordinary ASB interpretation.
- Exclude chrX from the ordinary analysis for female samples, or report a separate X-aware analysis.
- For cancer, supply measured per-variant RAF or matching gDNA through the tested BaalChIP interface; otherwise exclude altered segments and do not call the result copy-number-aware.
- Apply a matched genome blacklist and restrict tests to the called ChIP-seq peaks.

Keep the excluded set and reason codes so that biologically skewed loci are auditable rather than silently lost.

### 4. Run the selected model

For BaalChIP, use the shipped fail-closed single-group runner. Every assembly declaration must agree, the official five-column sample sheet and `ID/CHROM/POS/REF/ALT[/RAF]` het contract are validated, and the output directory must not already exist:

```bash
Rscript scripts/baalchip_workflow.R \
  --samples samples.tsv --hets tumor-hets.tsv --group TUMOR \
  --blacklist hg38-blacklist.bed --imprinted imprinted-hg38.bed \
  --sex female --assembly GRCh38 --samples-assembly GRCh38 \
  --hets-assembly GRCh38 --blacklist-assembly GRCh38 \
  --imprinted-assembly GRCh38 --correction raf \
  --samtools /opt/conda/bin/samtools --out results/TUMOR
```

Use `--correction gdna --gdna-bam matched-normal.bam` instead when gDNA should produce RAF. The runner removes any het-table RAF in that mode because BaalChIP otherwise gives it priority. It preserves pre-model exclusions with reason codes, selects the requested group from the report list, writes correction provenance and structured terminal state, and refuses an existing output directory. `--preflight-only` validates contracts without claiming model execution.

For RASQUAL, use `scripts/rasqual_cohort.py`: `prepare` validates text matrices and invokes the pinned source's official `txt2bin.R`; `run` requires a complete one-row-per-feature manifest, checks native-double binary sizes, invokes tabix and RASQUAL without a shell, rejects non-converged or malformed 25-column rows, and applies Benjamini-Hochberg correction to the rows actually emitted across the cohort run.

AlleleSeq is external-only in this skill. Do not execute a partial local sketch or substitute an unofficial `vcf2diploid` artifact. Follow the exact unavailable-toolchain boundary in [methods-and-commands.md](references/methods-and-commands.md), then hand off to a restored official environment.

### 5. Validate the result

At minimum:

- compare the genome-wide or tested-site REF allelic-ratio distribution with 0.5 after filtering;
- count input, remapped, retained, tested, and significant loci at each stage;
- inspect significant loci for blacklist, imprinting, chrX, and copy-number overlap;
- check replicate agreement and read depth at each reported site;
- verify REF/ALT orientation and chromosome naming in exported tables;
- distinguish measured ASB from sequence-model predictions such as chromBPNet or deepSEA.

A strong prediction without measured ASB can reflect insufficient coverage or model extrapolation; disagreement is not by itself proof that either result is wrong.

### 6. Report the analysis contract

Report the reference assembly, aligner, WASP or alternative bias strategy, genotype source, variant and peak filters, sample sex, chrX handling, imprinting resource, blacklist, copy-number source and handling, ASB method and version, model parameters, minimum depth, multiple-testing or posterior criterion, replicate strategy, and read/locus attrition.

Per-site output should retain chromosome, position, REF, ALT, allele counts, raw allelic ratio, corrected ratio when available, uncertainty or test statistic, significance call, depth, and exclusion/annotation flags.

## Interpretation boundaries

- ASB is a within-sample allelic measurement at heterozygous sites; a caQTL or bQTL is a population-level association.
- A significant imprinted or X-linked skew can be biologically real but is not automatically evidence of differential TF binding.
- Copy-number-driven allele dose can mimic binding imbalance.
- Personalized-genome alignment avoids reference bias but does not remove imprinting, X-inactivation, copy-number, low-depth, or peak-selection concerns.
- Deep-learning variant effects are predictions in a counterfactual sequence context; ASB measures the assayed sample and cell state.

## References

- Rozowsky J et al. 2011. *Molecular Systems Biology* 7:522 (AlleleSeq; <https://doi.org/10.1038/msb.2011.54>).
- van de Geijn B et al. 2015. *Nature Methods* 12:1061 (WASP; <https://doi.org/10.1038/nmeth.3582>).
- Kumasaka N et al. 2016. *Nature Genetics* 48:206 (RASQUAL; <https://doi.org/10.1038/ng.3366>).
- de Santiago I et al. 2017. *Genome Biology* 18:39 (BaalChIP; <https://doi.org/10.1186/s13059-017-1165-7>).
- Mayba O et al. 2014. *Genome Biology* 15:405 (MBASED).
- Chen J et al. 2016. *Nature Communications* 7:11101 (1000 Genomes ASB/ASE survey).

## Related skills

- chip-seq/peak-calling - call peaks upstream.
- chip-seq/chipseq-qc - evaluate ChIP-seq quality before ASB.
- chip-seq/chip-deep-learning - compare measured ASB with predicted variant effects.
- chip-seq/peak-annotation - annotate ASB variants to genes and cCREs.
- atac-seq/allele-specific-accessibility - parallel accessibility workflow.
- bio-causal-genomics-fine-mapping - use ASB as orthogonal fine-mapping evidence.
- bio-variant-annotation - annotate heterozygous variants.
- phasing-imputation/haplotype-phasing - prepare phased genotypes for AlleleSeq.
