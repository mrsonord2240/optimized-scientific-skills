# CRISPRi / CRISPRa TSS Targeting

## CRISPRi / CRISPRa TSS Targeting

**Goal:** Position guides relative to the empirical TSS for maximum knockdown (CRISPRi) or activation (CRISPRa).

**Approach:** Resolve TSS from FANTOM5 CAGE peaks (highest-ranked peak per gene; fall back to Ensembl/RefSeq if absent), define the modality-specific window, score candidate spacers in that window with Rule Set 2 plus the Horlbeck/Sanson CRISPRi/a-tailored rules, and select 5-6 guides per gene biased toward the window center.

Window functions and their rationale (Dolcetto vs Horlbeck v2 for CRISPRi; Calabrese vs Horlbeck v2 for CRISPRa; SAM vs SunTag) are in `scripts/tss_windows.py`:

```bash
python scripts/tss_windows.py --mode crispri --tss 1000000 --strand +   # prints start<TAB>end: 999950  1000300
python scripts/tss_windows.py --mode crispra --tss 1000000 --strand +   # 999850  999925
```

Caveat: the 75-bp Calabrese CRISPRa window routinely fails to hold a full 6-guide quota that passes the GC/poly-T filter. Report the actual count per gene rather than padding with out-of-window guides, or widen to Horlbeck v2 when the quota must be met.

**Critical nuance:** Cell-type-specific TSSs differ from the FANTOM5 consensus in ~15% of genes. For tissue-specific screens (e.g., neuron, hepatocyte), re-derive TSSs from a matched CAGE / GRO-seq / PRO-seq dataset before locking guide positions, or knockdown efficiency drops several-fold. The single most common cause of "weak" CRISPRi hits is mis-positioned guides against an alternative TSS: dCas9-KRAB knockdown is maximal within +/-100 bp of the actual Pol II loading site, and canonical Ensembl/RefSeq annotation can be off by 1-10 kb for genes with broad or non-canonical promoters. Symptom: "easy" essentials (RPS, RPL, EIF) show normal dropout but newer genes do not, and the library validates poorly against CEGv2.
