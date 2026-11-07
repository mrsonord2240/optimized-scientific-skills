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
