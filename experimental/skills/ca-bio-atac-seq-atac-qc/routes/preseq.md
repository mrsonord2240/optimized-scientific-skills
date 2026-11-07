# Library complexity extrapolation (preseq)

```bash
preseq c_curve -B -P -s 1e5 -o sample.ccurve.tsv sample.bam
preseq lc_extrap -B -P -e 200000000 -s 5000000 -o sample.lcextrap.tsv sample.bam
```

Shown, not scripted. Input must be a coordinate-sorted, duplicate-retained BAM (a deduplicated BAM has nothing to extrapolate). Paired-end ATAC needs `-P`; without it mates count as independent reads. For `c_curve` keep `-s` well below the fragment count (a step above depth returns only the `0 0` row; about a tenth of it works). preseq's default `-e` is 1e10 and default step 1M.

- If distinct fragments flatten before 100M reads, the library is bottlenecked: re-prep, because re-sequencing will not help.
- If the curve keeps climbing, deeper sequencing recovers more unique reads.

Done when the answer states re-prep or re-sequence with the curve as the reason.
