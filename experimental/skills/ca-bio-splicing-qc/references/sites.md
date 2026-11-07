## Splice Site Strength (MaxEntScan and SpliceAI)

**Goal:** Score donor and acceptor sites to flag weak / cryptic sites.

```python
import sys; sys.path.insert(0, 'scripts')
from splicing_qc import score_splice_sites

# 5'ss: 9 nt (3 exon + 6 intron). 3'ss: 23 nt (20 intron + 3 exon); the intron must end in AG.
donors = ['CAGGTAAGT', 'CAGATAAGT']
acceptors = ['TTTTTTTTTTTTTTCCTTAGGAG']          # 11.58; 'T'*20 + 'CAG' has no AG and scores -7.20
s5, s3 = score_splice_sites(donors, acceptors)   # one score per input; NaN for invalid input
print(s5, s3)                                     # [10.86, 2.68] [11.58]
```

`maxentpy.maxent.score5/score3` end the process (`SystemExit: Wrong length of fa!`) on a wrong length and raise `KeyError` on N/U; lower-case is accepted. The helper validates length and A/C/G/T first and returns NaN.

| Score | Interpretation |
|-------|----------------|
| 5'ss MaxEnt > 8 | Strong donor |
| 5'ss MaxEnt 5-8 | Moderate |
| 5'ss MaxEnt < 5 | Weak / cryptic |
| 3'ss MaxEnt > 8 / < 5 | Strong / weak acceptor |

MaxEntScan (Yeo & Burge 2004 *J Comput Biol*) defines the score; the cut-offs are conventions, supported on a chrX sample: 60.7% of 8,549 annotated GT donors scored >8 and 11.6% <5 (median 8.6), while 95.5% of 88 decoy GT donors scored <5 (AUC 0.97 donor, 0.96 acceptor), and all 399 non-GT annotated donors scored <5. Sampling and annotation release affect these descriptive figures; use them as rough context, not calibration. `score3` is comparatively slow (about 0.1 s/site in one checked environment), so scoring 10,000 acceptors can take roughly 20 minutes on one core; parallelize independent batches if needed.

**SpliceAI** predicts in-vivo usage from the full pre-mRNA context. Score a VCF (the assembly must match the FASTA):

```bash
spliceai -I variants.vcf -O variants.spliceai.vcf -R genome.fa -A grch38 -D 50 -M 0
# INFO: SpliceAI=ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL
```

Checked on a canonical donor G>A (PLCXD1, GRCh37): `DS_DG 0.90, DS_DL 1.00`. Published reference cut-offs on the maximum delta score are >=0.2 and <=0.1; 0.5 and 0.8 are the recommended and high-precision tiers of Jaganathan 2019 *Cell*. These are research prioritization thresholds here, not a clinical classification rule (Scope).

- **MaxEntScan** scores intrinsic sequence strength; **SpliceAI** scores contextual usage. High MaxEnt with low SpliceAI = intrinsically strong but contextually silenced; low MaxEnt with high SpliceAI = weak but contextually used (e.g. enhancer-driven). Report both; for variant impact see `splice-variant-prediction`.
