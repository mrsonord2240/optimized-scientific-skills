# Review context predictions
~~~bash
python epridict_prediction.py single --chromosome "$CHROM" --position_hg38 "$POSITION"
~~~
Set CHROM and POSITION from an existing input; run from the external ePRIDICT repository root. ePRIDICT is a version-dependent workflow requiring its model and chromatin tracks; the source has no bundled generic wrapper. It reports a model-specific locus score and percentile, not a per-guide prediction. Treat source context as conditional evidence rather than a causal conclusion; verify track integrity and assembly. The provided candidate example uses simulated scores rather than PRIDICT2 and cannot substantiate predictions.
Done when model, reference, input provenance and limitations are recorded or this analysis is deferred.
