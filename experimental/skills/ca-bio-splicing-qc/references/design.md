## Experimental Design Audit (Before Sequencing)

Targets are conventions, not hard limits.

| Decision | For splicing analysis | Rationale |
|----------|------------------------|-----------|
| **Library prep** | rRNA depletion (Ribo-Zero, RiboCop) | poly(A) selection loses pre-mRNA, nascent transcripts and detained introns; mandatory for IR analysis |
| **Read length** | PE 100-150 nt (PE 150 preferred) | Junction-spanning reads need >=8 nt overhang on each side |
| **Pairing** | Paired-end | Single-end loses fragment-level disambiguation of junctions |
| **Depth** | 50-100M reads/sample; >=100M for low-abundance events | 30M is the DGE-grade minimum; only a fraction of reads span junctions (24.5% on the chrX test data: 24,692 of 100,826 records) |
| **Strandedness** | Stranded library (dUTP / TruSeq stranded) | Distinguishes overlapping antisense; some tools double-count unstranded junctions |
| **Replicates** | n>=3 per condition | n=2 vs n=2 is poorly calibrated |
| **Long-intron genes (TTN, brain)** | Increase `--alignIntronMax` | Default 0 means the window-derived maximum, 2^16 x 9 = 589,824 nt (~590 kb; from `STAR --help`), which misses longer introns |

Annotation choice is in its own section below.
