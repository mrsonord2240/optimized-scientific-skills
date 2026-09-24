# Failure Modes

## Failure Modes

### Library skew from PCR bias during amplification

**Trigger:** Amplifying the cloned plasmid pool with too many PCR cycles (>20) or with high-GC-bias polymerase.
**Mechanism:** GC-extreme guides amplify nonlinearly; high-GC guides dominate, low-GC guides drop out.
**Symptom:** Gini >0.2 on plasmid pool; sgRNAs with GC <30% systematically depleted.
**Fix:** Cap PCR at 15 cycles; use Q5 or NEBNext Ultra II (low-bias); sequence at 500x post-amp to confirm Gini.

### Oligo-synthesis dropouts in low-complexity guides

**Trigger:** Chip-synthesis errors at homopolymer runs or guides starting with GGGG.
**Mechanism:** Synthesis chemistry has higher error rate at low-complexity regions; missing oligos cannot be cloned.
**Symptom:** Specific guides absent from plasmid pool despite no design-rule violation.
**Fix:** Re-design replacement guides; for production runs, request 2-3x synthesis depth so dropouts are buffered.

### Polyclonality from high MOI

**Trigger:** Infection at MOI >0.5 to "save cells."
**Mechanism:** Poisson math: at MOI 0.3, 26% of all cells are infected and 4% carry >=2 sgRNAs (14% of the infected fraction); at MOI 0.5, 39% are infected and 9% carry >=2.
**Symptom:** Hits include neutral genes that co-infect with true essentials.
**Fix:** MOI 0.3 strict; titer Cas9-positive cells specifically; re-check by qPCR of integration.
