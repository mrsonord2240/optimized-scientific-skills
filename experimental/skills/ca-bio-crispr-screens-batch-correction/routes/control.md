# Control
~~~bash
python scripts/ntc_normalize.py counts.tsv ntcs.txt normalized.tsv scale_factors.tsv
~~~
Supply a feature-indexed count matrix containing only sample columns and one control ID per line. The source recommends at least 500 controls and uses target median 1000. Zero or missing control medians cannot yield valid scaling. See [NTC notes](../references/ntc-anchored-normalization.md).
Done when the normalized matrix and scale factors are saved and control coverage checked.
