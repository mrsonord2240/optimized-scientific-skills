# Differential



```bash
# Compare control vs treatment nucleosome positions.
# The sample pair is the POSITIONAL argument (a:b means a minus b); -b is for background/input to
# subtract, and -c specifies a read-count to normalize to (an integer, NOT a control BAM path).
danpos dpos condition2.bam:condition1.bam \
    -o danpos_diff/ \
    --paired 1 \
    --smooth_width 80
```

The differential table is `condition2-condition1.positions.integrative.xls`. Its shift column is `treat2control_dis` (signed bp; `diff_smt_loca` is a coordinate) and the FDR that tracks position shifts is `point_diff_FDR`. `smt_diff_FDR` tests summit height and removed every planted shift in testing. The table has no event-type label column; classify gained, lost, shifted, and fuzziness-changed nucleosomes from the value and FDR columns yourself. ENCODE has no official threshold; require |shift| >= 30 bp and `point_diff_FDR` < 0.05:

```bash
awk -F'	' 'NR==1{for(i=1;i<=NF;i++)c[$i]=i; print; next}
     $c["treat2control_dis"]!="NA" && $c["point_diff_FDR"]!="NA" &&
     ($c["treat2control_dis"]>=30 || $c["treat2control_dis"]<=-30) && $c["point_diff_FDR"]<0.05'     danpos_diff/condition2-condition1.positions.integrative.xls > shifted.tsv
```

Full ATAC-tuned DANPOS3 recipe:

```bash
# --width 145: summit-scan window (DANPOS -jw/--width; default 40)
# --smooth_width 80: smoothing kernel width (DANPOS -z; default 20, widened for ATAC)
# -jd 145: min distance between adjacent nuc calls (single-dash short flag)
# --pheight 1e-5: occupancy P-value cutoff (DANPOS -p; dpos default 0). -q/--height is the separate density cutoff
# --frsz 200: fragment size used (mono-nuc)
danpos dpos sample.bam \
    --paired 1 \
    --width 145 \
    --smooth_width 80 \
    -jd 145 \
    --pheight 1e-5 \
    --frsz 200 \
    --out danpos_out/
```

Verify exact flags with `danpos dpos --help`. The tested route is the conda `danpos3` package (3.2.4), which installs an executable named `danpos` and prints a 3.1.1 banner; a GitHub checkout (github.com/sklasfeld/DANPOS3) runs as `python danpos.py` with the same flags. The bioconda `danpos` package is DANPOS2.

Adapted from DANPOS3 docs for ATAC; `--smooth_width 80` widens the smoothing kernel to match ATAC's sharper signal versus MNase's broader cleavage. `-jd 145` (alternative `--distance 145`) enforces nucleosome spacing >= 145 bp (one nucleosome footprint).
