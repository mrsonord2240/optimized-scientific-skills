## AP-MS / Affinity-Enrichment Scoring

**Goal:** Rank prey in a pulldown by enrichment over NEGATIVE-CONTROL IPs, not by abundance or by ratio to the input lysate.

**Approach:** A pulldown is deliberately non-representative, so no data-internal normalization (median, sample-loading, IRS) applies, and the input lysate is not a control -- sticky background (ribosome, chaperones, tubulin, keratin) binds the beads in the pulldown and is diluted in the input, so "top N over input" returns background. The control IP is the only thing that separates a bead binder from an interactor. Require reproducible detection across bait replicates and enrichment over control, and floor absent controls at the run's detection limit so bait-only prey score finitely instead of `+Inf`. The default `min_bait_reps` is every bait replicate (reproducibility first), which trades sensitivity for specificity: at 3 replicates a true interactor missing from one by stochastic dropout is dropped (an enrichment of 5.4 was, on the audit fixture). `min_bait_reps=2` is the usual compromise; run both and report how many prey the looser setting adds so the choice is visible.

```bash
python scripts/apms_score.py ip.csv --bait Bait1,Bait2,Bait3 --ctrl Ctrl1,Ctrl2,Ctrl3 [--min-bait-reps 2] --out scores.csv
```
Or `from apms_score import score_vs_control_ips`. Input is prey x replicate RAW intensities or spectral counts (0/NaN = not detected); prey never seen in any bait IP get a NaN enrichment and are never called. Output columns: `n_bait`, `n_ctrl`, `n_ctrl_runs_used`, `log2_enrichment`, `worst_case_log2` (weakest bait replicate vs best control), `interactor`. A control IP with no data is printed and excluded; all controls empty raises.

Seeded matrix (sticky binders excluded, dead control announced, `min_bait_reps` dropout): `examples/lfq_normalization.py`.

This fold-change/presence score is the honest ceiling for ONE bait with a few controls. For a
probability rather than a cutoff use SAINTexpress (spectral counts, `interaction`/`prey`/`bait`
files -> AvgP, report BFDR <= 0.01-0.05) or CompPASS WD scores across a bait MATRIX; both need
several independent baits, or the CRAPome as an external control set, before their statistics mean
anything. Feed them counts or raw intensities -- never a median/SL/IRS-normalized matrix.
