# Combat
~~~bash
python scripts/combat_correct.py counts.tsv metadata.tsv corrected.tsv uncorrected.tsv --batch-col batch --condition-col condition
~~~
Check sample matching, nonnegative counts, no NaNs, at least three samples per batch, and a non-confounded design. The source script uses log2(counts+1), models condition and leaves features with insufficient residual variance uncorrected. Review the uncorrected list and compare biological and batch structure. See [ComBat notes](../references/combat.md).
Done when the corrected matrix and uncorrected-feature list are saved and reviewed.
