# Inspect SILAC data
~~~bash
python scripts/silac_checks.py pilot_evidence.txt > pilot_checks.json
~~~
The source expects a heavy-only pilot with Intensity H/L and, for Arg-to-Pro, Sequence and Ratio H/L. Labeling >=95% is the source convention; report missing sequence and slope limitations. To retain on/off cases separately from finite ratios:
~~~bash
python scripts/silac_ratio.py peptide_table.tsv silac_ratios.tsv --heavy "Intensity H" --light "Intensity L"
~~~
Both channels positive yield log2 H/L; H-only and L-only are presence flags, not finite fold changes. See [SILAC notes](../references/silac.md).
Done when pilot checks, ratio presence and limitations are reviewed.
