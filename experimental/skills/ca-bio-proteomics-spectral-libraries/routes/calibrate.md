# Calibrate retention times
~~~bash
python scripts/library_tools.py calibrate anchors.tsv library.tsv calibrated.tsv
~~~
Anchor TSV needs irt and observed_rt; library needs iRT. The source linear fit requires R-squared >=0.95 and reports slope and intercept; nonlinearity or inadequate anchors require separate evaluation, potentially LOWESS. Confirm observed RT units, anchor identity, model and acquisition. The output adds calibrated_rt rather than assuming a destination engine's column convention.
Done when fit diagnostics, units, output and alignment assumptions are reviewed.
