# IDR-passed union or quick merge

IDR-filtered union (reproducibility-aware, ENCODE style): run IDR per true replicate pair, then

```bash
cat pair1.idr.bed pair2.idr.bed | sort -k1,1 -k2,2n | bedtools merge -i - > idr_union.bed
```

- Use true-replicate IDR at 0.05 only, or pseudoreplicate IDR at 0.10 only. Do not mix; pseudoreplicate peaks are looser.
- The IDR set has far fewer peaks than overlap counting; use it for high-confidence reporting.

Exploratory merge (no statistics, widths inflate):

```bash
cat rep*.narrowPeak | sort -k1,1 -k2,2n | bedtools merge -i - > merged.bed
```

Intersection across all replicates (`bedtools multiinter`, keep rows where the count column equals the number of samples) is high stringency but drops condition-specific peaks.

Both commands are shown, not scripted. The result is variable width: not for differential counting.

Done when the answer states the union rule and the IDR threshold used; next `routes/count.md`.
