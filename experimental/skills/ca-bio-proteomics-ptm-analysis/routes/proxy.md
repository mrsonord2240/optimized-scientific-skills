# No-global proxy analysis
~~~bash
python scripts/run_proxy_checked.py --timeout 1800 --grace 0.15 --receipt "$RECEIPT" --stage-dir "$STAGE" --publish-dir "$FINAL_DIR" -- "$RSCRIPT" scripts/msstatsptm_labelfree.R "dir=$DATA" "out=$STAGE" use_unmod=TRUE
~~~
Use fresh stage and final directory and set all path variables. Paired-global files must be absent; the worker validates sample annotation, unambiguous peptide/protein assignment, unmodified-peptide coverage and replicates. Accept the bundle only when its internal manifest and marker-last proxy_checked_complete.json both validate, with clean worker and no lingering descendants. Outputs are proxy-adjusted candidates, not protein-adjusted regulation calls; retain raw PTM results.
Done when receipt, marker, hashes and non-regulatory labels are verified.
