# Prime editor quantification (CRISPResso2)

Moved from SKILL.md. Read when quantifying pegRNA-templated edits.

## Prime Editor Quantification

**Goal:** Quantify pegRNA-templated edits versus indel byproducts and partial edits.

**Approach:** Provide spacer, extension (PBS + RTT), and scaffold sequences; CRISPResso identifies reads matching the intended edit.

```bash
CRISPResso \
    --fastq_r1 pe_sample.fastq.gz \
    --amplicon_seq <amplicon_seq> \
    --guide_seq <20nt_protospacer> \
    --prime_editing_pegRNA_spacer_seq <20nt_protospacer> \
    --prime_editing_pegRNA_extension_seq <RTT+PBS_sequence> \
    --prime_editing_pegRNA_scaffold_seq <scaffold_sequence> \
    --output_folder pe_results \
    --name pe_sample

# Output adds:
#   Prime-editing outcomes are extra amplicon rows (Reference / Prime-edited / Scaffold-incorporated)
#   inside CRISPResso_quantification_of_editing_frequency.txt
```

**Reading prime-editor output:**

| Metric | Interpretation |
|--------|----------------|
| Intended edit % | The pegRNA-encoded edit was correctly installed |
| Scaffold incorporation % | Reverse transcription read into scaffold instead of stopping at edit; failure mode |
| Indel % | Nick-only editing without templated repair; common at low-PE-activity sites |
| Unmodified % | Read matches the reference exactly |

Acceptance levels are in Quantitative Thresholds. See [[prime-editing-screens]] for pegRNA design rules.

### Prime editor sample with high scaffold incorporation

**Trigger:** RTT is too short relative to PBS, or pegRNA stops short.
**Mechanism:** Reverse transcriptase reads past the edit into scaffold sequence; product is detectable but undesired.
**Symptom:** Scaffold incorporation >5%; intended edit efficiency lower than expected.
**Fix:** Re-design pegRNA with longer RTT; verify with PRIDICT2 (see [[prime-editing-screens]]).
