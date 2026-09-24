# CBE / ABE quantification (CRISPResso2)

Moved from SKILL.md. Read when quantifying a cytosine or adenine base-editing experiment.

## Base Editor Quantification

**Goal:** Distinguish target base conversion from bystander edits and indel byproducts.

**Approach:** Run CRISPResso with `--base_editor_output` flag and specify the conversion direction; widen the quantification window to cover the editing window.

```bash
# Cytosine Base Editor (CBE): C->T conversion
CRISPResso \
    --fastq_r1 cbe_sample.fastq.gz \
    --amplicon_seq <amplicon_seq> \
    --guide_seq <20nt_protospacer> \
    --base_editor_output \
    --conversion_nuc_from C \
    --conversion_nuc_to T \
    --quantification_window_size 10 \
    --quantification_window_center -10 \
    --output_folder cbe_results \
    --name cbe_sample

# Adenine Base Editor (ABE): A->G conversion
CRISPResso \
    --fastq_r1 abe_sample.fastq.gz \
    --amplicon_seq <amplicon_seq> \
    --guide_seq <20nt_protospacer> \
    --base_editor_output \
    --conversion_nuc_from A \
    --conversion_nuc_to G \
    --quantification_window_size 10 \
    --quantification_window_center -10 \
    --output_folder abe_results \
    --name abe_sample
```

**Reading the output:**

| Metric | Where | Interpretation |
|--------|-------|----------------|
| Target editing % | `Quantification_window_nucleotide_percentage_table.txt`, target C/A row | Primary endpoint |
| Bystander editing % | Same table, other C/A positions in window | Off-target byproduct in window |
| Indel rate | `CRISPResso_quantification_of_editing_frequency.txt` | Cas9-like cut artifacts (limit in Quantitative Thresholds) |
| Substitution-vs-indel ratio | Derived | Distinguishes clean BE from cut-mediated mutagenesis (cutoffs in Quantitative Thresholds) |

**Critical:** Bystander editing is intrinsic to base editors (the deaminase acts across a 5-nt window); it is not noise. Report bystander rates alongside target rates. See [[base-editing-analysis]] for variant-call implications.

### Bystander C/A editing inflates "editing efficiency"

**Trigger:** Base-editor sample with wide quantification window; bystander Cs at adjacent positions counted as edits.
**Mechanism:** A wide window (`--quantification_window_size 10`; the CRISPResso default is 1) includes all positions in the editing window; bystander edits are real but distinct from target edit.
**Symptom:** Editing efficiency 80%+ but target SNV is 30%; bystander rate is 50%.
**Fix:** Always read the per-position table (`Quantification_window_nucleotide_percentage_table.txt`), not just the aggregate. Report target and bystander rates separately. See [[base-editing-analysis]].
