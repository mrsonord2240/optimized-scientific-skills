# Paired-global protein adjustment
~~~bash
python scripts/run_checked.py --timeout 1800 --grace 0.15 --receipt "$RECEIPT" --stage-dir "$STAGE" --publish-source adjusted_sites.csv --publish-dest "$FINAL" --csv adjusted_sites.csv --require-columns Protein,Label,log2FC,adj.pvalue_lfc --min-rows 1 -- "$RSCRIPT" scripts/msstatsptm_labelfree.R "dir=$DATA" "out=$STAGE" use_unmod=FALSE
~~~
Set DATA, RSCRIPT, STAGE, RECEIPT and FINAL to appropriate paths; stage and final must be fresh. Supply matched enriched/global inputs and annotations. Checked runner requires clean worker exit, no owned descendants and valid declared CSV before no-clobber publication; bare R execution is diagnostic-only. Source uses class-I 0.75, Treatment vs Control and lfc=1 in its TREAT-style test. TMT uses the separate checked command and full output schema in [TMT](../references/tmt-isobaric.md), including a reference channel per plex.
Done when receipt, output schema, contrast, paired design and assumptions are reviewed.

TMT worker: `scripts/msstatsptm_tmt.R`. The optional one-process Windows teardown guard is `scripts/cli_windows_cleanup_guard.Rprofile`; use only as specified by source runtime guidance.
