---
name: bio-crispr-screens-crispresso-editing
category: Data Analysis
description: Quantifies CRISPR editing outcomes with CRISPResso2 (Clement 2019 Nat Biotechnol) across Cas9-nuclease (indels, HDR), CBE and ABE base editors (target conversion + bystander), and prime editor (pegRNA-templated) modes. Covers single-amplicon (CRISPResso), multi-sample batch (CRISPRessoBatch), pooled-amplicon (CRISPRessoPooled), WGS off-target (CRISPRessoWGS), and sample-comparison (CRISPRessoCompare) workflows; quantification-window math that controls what is called edited; substitution-vs-indel diagnostic to distinguish BE from Cas9 contamination; MMEJ deletion pattern interpretation; allele-frequency tables; and failure modes from amplicon misalignment or contamination. Use when quantifying editing from amplicon sequencing, choosing CRISPResso mode by design, distinguishing intended edits from bystanders and indel byproducts, debugging low-alignment runs, or generating publication-grade editing reports.
tool_type: cli
primary_tool: CRISPResso2
license: MIT
author: GPTomics
---

## Version Compatibility

Reference examples checked on CRISPResso2 2.3.4 (pinellolab/crispresso2 Docker image), pandas 2.2+, numpy 1.26+, matplotlib 3.8+.

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `CRISPResso --version`; `CRISPRessoBatch --help`; `CRISPRessoPooled --help`; `CRISPRessoWGS --help`; `CRISPRessoCompare --help`
- Python: `from CRISPResso2 import ...`

If code throws ImportError, AttributeError, or TypeError, introspect the installed package and adapt the example to match the actual API rather than retrying.

Install: `conda install -c bioconda crispresso2` (not on PyPI; bioconda has no win-64 build, so on Windows use the
`pinellolab/crispresso2` Docker image). Plotting/parsing helpers: `pip install pandas matplotlib seaborn`.

Required inputs:
- FASTQ files (paired-end recommended, single-end accepted)
- Amplicon reference sequence, primer-trimmed (do NOT include primer regions)
- Guide protospacer sequence (20 nt; PAM not included)
- For HDR: expected edited amplicon sequence
- For BE: target nucleotide conversion direction (C->T for CBE; A->G for ABE)
- For PE: pegRNA spacer + extension (RTT+PBS) + scaffold sequences

## CRISPResso2 Editing Quantification

**"Quantify CRISPR editing from my amplicon sequencing"** -> Align amplicon reads against the reference, classify each read as unmodified / NHEJ / HDR / base-edited / prime-edited within the quantification window, and report per-edit-type frequencies, indel size distributions, allele-frequency tables, and substitution-position profiles.

- CLI: `CRISPResso` -- single amplicon, single sample
- CLI: `CRISPRessoBatch` -- multi-sample with per-sample parameters
- CLI: `CRISPRessoPooled` -- multi-amplicon pooled amplicon sequencing
- CLI: `CRISPRessoWGS` -- off-target quantification from whole-genome BAM
- CLI: `CRISPRessoCompare` -- pairwise outcome comparison (e.g., treated vs untreated)

## Mode Decision Tree

| Experimental design | Mode | Key parameters |
|---------------------|------|----------------|
| Single amplicon, single sample (e.g. pilot edit validation) | `CRISPResso` | `--amplicon_seq`, `--guide_seq` |
| Same amplicon, many samples (e.g. timecourse, dose response) | `CRISPRessoBatch` | `--batch_settings` table; see `references/batch-pooled-wgs.md` |
| Many amplicons, pooled in one library (e.g. arrayed validation pool) | `CRISPRessoPooled` | `--amplicons_file`; read `references/batch-pooled-wgs.md` first (default `--min_reads_to_use_region` silently yields all-`NA` rows) |
| Off-target survey from whole-genome BAM | `CRISPRessoWGS` | `--bam_file`, `--reference_file`, `--region_file`; read `references/batch-pooled-wgs.md` (low-read regions come back `NA`) |
| Comparing two CRISPResso runs (e.g. condition A vs B) | `CRISPRessoCompare` | two positional output folders |
| HDR / knock-in validation | `CRISPResso` with `--expected_hdr_amplicon_seq` | Same as base CRISPResso |
| Cytosine base editor (C->T) | `CRISPResso --base_editor_output` | `--conversion_nuc_from C --conversion_nuc_to T`; see `references/base-editor.md` |
| Adenine base editor (A->G) | `CRISPResso --base_editor_output` | `--conversion_nuc_from A --conversion_nuc_to G`; see `references/base-editor.md` |
| Prime editor (templated edit) | `CRISPResso` with pegRNA parameters | `--prime_editing_pegRNA_spacer_seq`, `--prime_editing_pegRNA_extension_seq`, `--prime_editing_pegRNA_scaffold_seq`; see `references/prime-editor.md` |

**Fails when:**
- Pooled-amplicon mode applied to amplicons that share primer sequences -- reads get assigned to whichever amplicon comes first. Design primers with >=3-bp distinguishing regions or use unique molecular identifiers.
- Base editor mode without specifying `--conversion_nuc_from`/`--conversion_nuc_to` -- defaults assume CBE (C->T); ABE runs will misclassify.
- Prime editor mode without `--prime_editing_pegRNA_extension_seq` -- the RTT template is missing, no edit is detectable.

## Reference Files

| File | Read when |
|------|-----------|
| `references/base-editor.md` | Quantifying CBE/ABE editing: worked commands, target vs bystander reading table, bystander-inflation failure mode |
| `references/prime-editor.md` | Quantifying pegRNA-templated edits: worked command, intended-edit / scaffold reading table, high-scaffold failure mode |
| `references/batch-pooled-wgs.md` | Running CRISPRessoBatch, CRISPRessoPooled or CRISPRessoWGS: settings-file formats, output layout, the `NA`-row thresholds |

## The Quantification Window

**Why this matters for postdoc-level use:** CRISPResso classifies reads as "edited" or "unmodified" based on whether *modifications fall inside the quantification window* (not the whole amplicon). The window is centered on the predicted cut site (Cas9: 3 bp upstream of PAM; Cas12a: 18 bp downstream of PAM) with a default size of 1. `--quantification_window_size N` extends N bp on EACH side, so the window is 2N bp wide.

```bash
# Default Cas9 setup
--quantification_window_size 1                  # 1-bp window at cut site
--quantification_window_center -3               # 3 bp upstream of PAM

# Base editor: widen window to cover editing positions 4-8
--quantification_window_size 10                 # 10 bp each side (20 bp total)
--quantification_window_center -10              # center on the editing window
```

**Consequences of mis-sized window:**
- Too narrow: misses edits at HDR positions or far bystanders; underestimates editing
- Too wide: includes random sequencing errors; inflates editing rate
- Wrong center: edits at correct position are scored as outside the window

For base editing screens a widened window is conventional; CRISPResso2's own base-editor guidance uses `--quantification_window_center -17` with a window sized to span the editing positions. For prime editing with multi-base templated edits, widen to encompass the entire edit region.

## Single-Amplicon Cas9 Editing

**Goal:** Quantify indel frequencies and HDR efficiency from a single target site.

**Approach:** Align FASTQ reads to the reference and (optional) expected-HDR amplicon, classify each read, and report aggregated statistics.

```bash
# --expected_hdr_amplicon_seq is optional (HDR only); --min_average_read_quality is a Phred filter (see note below)
CRISPResso \
    --fastq_r1 sample_R1.fastq.gz \
    --fastq_r2 sample_R2.fastq.gz \
    --amplicon_seq <amplicon_sequence_ref_genome> \
    --guide_seq <20nt_protospacer_no_PAM> \
    --expected_hdr_amplicon_seq <edited_amplicon_for_HDR> \
    --quantification_window_size 1 \
    --quantification_window_center -3 \
    --min_average_read_quality 30 \
    --output_folder sample_results \
    --name sample_id

# Outputs:
#   sample_results/CRISPResso_on_<name>/CRISPResso_mapping_statistics.txt
#   sample_results/CRISPResso_on_<name>/CRISPResso_quantification_of_editing_frequency.txt
#   sample_results/CRISPResso_on_<name>/Alleles_frequency_table.zip
#   sample_results/CRISPResso_on_<name>/3a.Indel_size_distribution.pdf
#   sample_results/CRISPResso_on_<name>/4b.Insertion_deletion_substitution_locations.pdf
#   (PDF by default; add --save_also_png for PNG)
```

**`--min_average_read_quality` changes the reported editing percentage** -- it is a real filter, not a
neutral default. On one real test amplicon it dropped aligned reads 235->221 and shifted Modified% from
26.38% to 24.89%. Always report the retained-read fraction (`READS ALIGNED` / `READS IN INPUTS` after
filtering, from `CRISPResso_mapping_statistics.txt`) alongside the editing percentage so a reader can see
how much filtering happened.

**Key outputs:**

| File | Content |
|------|---------|
| `CRISPResso_mapping_statistics.txt` | Tab-separated, one data row: READS IN INPUTS, READS AFTER PREPROCESSING, READS ALIGNED, N_COMPUTED_ALN, ... (no percentage columns) |
| `CRISPResso_quantification_of_editing_frequency.txt` | % unmodified, % NHEJ, % HDR (if expected), per-edit-class breakdown |
| `Alleles_frequency_table.zip` | Per-allele sequences and frequencies (allele-level resolution) |
| `Nucleotide_percentage_table.txt` | Per-position A/C/G/T/- frequencies (substitutions + deletions) |
| `Quantification_window_nucleotide_percentage_table.txt` | Same, restricted to quantification window (base-editor analysis) |


## Parse Output in Python

**Goal:** Pull editing metrics into downstream analysis or reports.

**Approach:** Read the tab-separated quantification files and the JSON metadata.

```bash
python scripts/parse_crispresso.py sample_results/CRISPResso_on_sample_id   # metrics as JSON
python scripts/parse_crispresso.py --selftest                               # parser self-check
```

From Python, `from parse_crispresso import parse_crispresso` (with `scripts/` on `sys.path`) returns
`reads_in_input`, `reads_aligned`, `mapping_pct` (computed as `READS ALIGNED / READS IN INPUTS`, since the
mapping-statistics file has no percentage column), `editing_quant` and, when present, `info`
(`CRISPResso2_info.json`). Checked on the FANC.Cas9 test amplicon: 250 / 235 / 94.0%, Modified% 26.38.

## Failure Modes

### Low alignment rate (<50%)

**Trigger:** Wrong amplicon sequence (off by one nt, wrong strand, primer-trimmed vs untrimmed).
**Mechanism:** CRISPResso fails to align reads beyond the amplicon edges; discards as unmappable.
**Symptom:** `READS ALIGNED` / `READS IN INPUTS` (from `CRISPResso_mapping_statistics.txt`, see `scripts/parse_crispresso.py`) <50%; per-position coverage drops at amplicon edges.
**Fix:** Re-derive amplicon from genome at primer-trimmed boundaries; verify strand orientation; check that primers are NOT included in `--amplicon_seq`.

### Total alignment failure (wrong locus / zero output)

**Trigger:** Amplicon sequence is from the wrong locus entirely (e.g. wrong gene), not just off-by-a-few-nt.
**Mechanism:** No reads align at all; CRISPResso exits with a hard error instead of writing a graded result.
**Symptom:** `CRITICAL: Alignment error, please check your input. / ERROR: No alignments were found`, exit code 1, no output folder written -- there is no percentage file to inspect in this case.
**Fix:** Confirm the amplicon sequence actually corresponds to the intended target locus (BLAT/BLAST it against the reference genome) before re-checking primer trimming or strand.

### High substitution rate but low indel (Cas9 sample)

**Trigger:** Sample contamination with adjacent amplicon, primer-dimer, or sequencing error inflation.
**Mechanism:** Random substitutions inflate the per-position substitution rate without true indels.
**Symptom:** Substitutions >2% at base positions outside the cut site; alignment metrics look fine.
**Fix:** Increase `--min_average_read_quality` to 30+; filter contaminating amplicons; check primer-dimer in `CRISPResso_RUNNING_LOG.txt`.

### MMEJ deletion misclassified as NHEJ

**Trigger:** Deletions with microhomology at junction; CRISPResso reports them as indels but doesn't distinguish MMEJ.
**Mechanism:** MMEJ creates predictable deletions using flanking microhomologies; biologically distinct from random NHEJ.
**Symptom:** Recurring same-size deletions in allele table (e.g., -7 bp deletion in 30% of reads).
**Fix:** Examine `Alleles_frequency_table` for over-represented allele patterns; flag MMEJ-mediated deletions for interpretation (these may be inferred from indel hotspots).

## Quantitative Thresholds

| Threshold | Value | Source / Rationale |
|-----------|-------|--------------------|
| Cas9 editing efficiency (functional KO) | >70% indels | Field convention; below this, KO is incomplete |
| Indel rate (clean base editor) | <5% | Field convention; >5% = unwanted cut activity (usually Cas9 / nCas9 expression mismatch) |
| Substitution-vs-indel ratio (BE) | >10 = clean BE; <3 = cut-mediated (Cas9-like) mutagenesis, use Cas9-like analysis | Diagnostic for BE purity |
| Target conversion (CBE) | >30% | Variable by target; below this, screen power is poor |
| Target conversion (ABE) | >30% | ABE typically lower per-base than CBE |
| Bystander rate (BE) | <10% acceptable; <5% ideal | Application-dependent; for variant function studies, must be controlled |
| Intended-edit % (prime editor) | >5% per-edit | Field convention; >20% at favorable sites, can be 50%+ |
| Scaffold incorporation (PE) | <2% | High-quality pegRNA design |
| Alignment rate | >85% | Below this, amplicon design or contamination issue |
| Minimum read quality | Phred 30 | Q30 Illumina base-call-accuracy standard |
| Read depth per sample | 1,000+ | Reliable allele table; higher for low-frequency variants |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| All reads "modified" | Misaligned reference | Check amplicon strand; reverse-complement test |
| BE shows mostly indels | Cas9 contamination or wrong protein | Re-derive cell line origin; check Cas9 vs nCas9-BE3 |
| Inconsistent batch results | Different amplicon_seq per sample | Use CRISPRessoBatch with consistent amplicon |
| Out-of-window edits ignored | Window too narrow | Increase `--quantification_window_size` |
| Allele frequency dominated by 1 read | Low input / clonal | Verify input cell count; rerun if singleton |

## References

- Clement K et al. 2019. *Nat Biotechnol* 37:224. CRISPResso2 algorithm and modes.
- Pinello L et al. 2016. *Nat Biotechnol* 34:695. Original CRISPResso.
- Anzalone AV et al. 2019. *Nature* 576:149. Prime editing (PE-1/PE-2/PE-3).
- Komor AC et al. 2016. *Nature* 533:420. Base editing (BE3).
- Findlay GM et al. 2018. *Nature* 562:217. Saturation genome editing.

## Related Skills

- crispr-screens/base-editing-analysis - Variant-function analysis using CRISPResso2 BE output
- crispr-screens/prime-editing-screens - PRIDICT2 pegRNA design + PE-tiling
- crispr-screens/library-design - sgRNA / pegRNA design for editing screens
- crispr-screens/screen-qc - Editing-efficiency QC for variant interpretation
- variant-calling/variant-annotation - Annotate detected variants downstream
- read-alignment/bwa-alignment - For WGS off-target alignment input
