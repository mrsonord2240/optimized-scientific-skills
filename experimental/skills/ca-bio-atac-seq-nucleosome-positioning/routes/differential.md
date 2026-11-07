# Differential nucleosome positioning

```bash
danpos dpos condition2.bam:condition1.bam -o danpos_diff/ --paired 1 --smooth_width 80
```

The positional pair order is condition2 minus condition1; `-b` is background, not control BAM. Filter signed `treat2control_dis` by absolute shift >=30 bp and `point_diff_FDR` <0.05, as in `references/differential.md`; `smt_diff_FDR` tests height instead. Verify installed DANPOS3 flags and input format. Compare biological replicates and report assumptions and thresholds.

Done when the differential table and filtered shifts are checked and the contrast is stated.
