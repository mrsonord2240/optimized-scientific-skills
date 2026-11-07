# Background (not needed to run an analysis)

The 10.4 bp helical periodicity on the nucleosome-free peak reflects the helical pitch of B-form DNA, with Tn5 preferring outward-facing minor grooves on nucleosomal DNA. Its presence is a positive indicator, not a requirement.

Why the two TSS scores differ: the ENCODE-style score divides a 100 bp centre window by 100 bp flanks at +/-1900 to 2000 bp, per-base; ATACseqQC `TSSEscore` sums TSS +/-100 bp over +/-1000 bp flanking windows, so ratios are larger (typically 2-3x). deepTools `plotProfile` gives a profile with no standard numeric score and is for visualization only.

Why mitochondrial reads dominate: mtDNA is naked of histones, so Tn5 cuts it hyperactively unless nuclei are isolated.

Why PCR bottlenecking shows as low NRF: each PCR cycle doubles starting fragments; with under 5000 starting cells and more than 12 cycles, distinct fragments saturate and reads pile at identical positions. FRiP tracks TSS enrichment because both measure how concentrated the signal is.

Spearman correlation of 10 kb bins tracks coverage similarity; 0.9 is a working target for biological replicates.

Spike-in: normalization choice materially changes differential accessibility when a global shift is expected, which motivates an external chromatin reference.

The source Skill had a separate usage guide of prompts and a threshold quick-reference duplicating `routes/thresholds.md`; they were dropped.
