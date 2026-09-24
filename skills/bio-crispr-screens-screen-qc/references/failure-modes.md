# Failure Modes

### High Gini in plasmid pool despite passing all design rules

**Trigger:** Library was cloned and amplified through too many PCR cycles (>20) or used a high-GC-bias polymerase.
**Mechanism:** Each PCR cycle compounds GC bias by ~5%; high-GC and low-GC guides become non-linear functions of starting abundance.
**Symptom:** Gini >0.15 in plasmid, GC-content stratification of dropout.
**Fix:** Cap PCR at 15 cycles for amplification; use Q5 / NEBNext Ultra II / KAPA HiFi (low-bias); re-sequence post-amp; if still bad, re-clone from glycerol stock.

### Falling PR-AUC across timepoints despite stable Gini

**Trigger:** Cas9 was not selected for before screen start; Cas9-negative cells in the pool dilute essentiality signal.
**Mechanism:** Each Cas9-negative cell carries a sgRNA but no editing; its sgRNA persists despite biological essentiality of the target.
**Symptom:** PR-AUC declines from 0.7 at week 1 to 0.4 at week 3; Gini and Pearson both pass.
**Fix:** Always select Cas9-positive cells (FACS or blast) before infection. For a salvage of an already-run screen, model Cas9-expression heterogeneity as a noise floor and accept reduced sensitivity.

### Apparent essentiality of amplified loci

**Trigger:** Cancer cell line with focal amplification (ERBB2 in SK-BR-3, MYC in colorectal, FGFR1 in head and neck).
**Mechanism:** Aguirre 2016 / Munoz 2016: many simultaneous Cas9 cuts trigger a DNA-damage response and G2 arrest; sgRNAs at amplified loci appear depleted independently of target essentiality.
**Symptom:** Hits include genes within known amplicons; sgRNAs with more genome-wide cut sites are more depleted.
**Fix:** Apply CRISPRcleanR pre-hoc or use Chronos/CERES with matched CN profile (see [[copy-number-correction]]). Always required for cancer-cell-line screens, not optional.

### Outlier replicate dragging Pearson down

**Trigger:** One technical replicate had a library-prep failure (low input, PCR jackpot, sequencing-lane swap).
**Mechanism:** Outlier sample has different total reads or different per-sgRNA distribution but passes individual sample QC.
**Symptom:** Pearson between replicates 0.85-0.90 with one pair as outlier; condition-level means look fine.
**Fix:** Drop the outlier replicate; re-derive Pearson on the remaining pair. If only two replicates and one is outlier, the condition lacks replication and must be re-run.

### Low Day-0 coverage from high MOI

**Trigger:** Infection at MOI >0.5.
**Mechanism:** Multiple sgRNAs per cell (Poisson table under MOI Verification); the "single-perturbation" assumption underlying every analysis method is violated.
**Symptom:** Apparent gene-gene interactions in single-gene screens; gene-level z-scores noisy; Pearson lower than expected for high-quality counts.
**Fix:** Re-titrate, re-infect at MOI 0.3, re-run screen (no analytical correction, see MOI Verification decision rule).

### CRISPRi/a screen with no signal on validated essentials

**Trigger:** Library targets wrong TSS (Ensembl canonical vs FANTOM5 highest-rank).
**Mechanism:** dCas9-KRAB knockdown is maximal within ±100 bp of the actual Pol II loading site; canonical annotation can be off 1-10 kb.
**Symptom:** RPS/RPL/EIF families dropping out as expected (these have clean canonical TSSs) but downstream genes failing; PR-AUC on broader CEGv2 panel drops.
**Fix:** Re-design library against FANTOM5 highest-CAGE-peak TSS (Sanson 2018); for tissue-specific lines, use matched CAGE / GRO-seq.
