# Search with pyOpenMS
~~~bash
python scripts/pyopenms_search.py sample.mzML target_decoy.fasta search_results.idXML --precursor-ppm 10 --fragment-da 0.02 --missed-cleavages 2
~~~
Use a centroided mzML and matched target-plus-decoy FASTA; defaults include fixed carbamidomethyl C and variable oxidation M. Next annotate and filter:
~~~bash
python scripts/pyopenms_fdr.py search_results.idXML target_decoy.fasta psms_1pct.idXML --decoy-string DECOY_ --fdr 0.01
~~~
Ensure the prefix matches and both targets and decoys are recognized. Check pyOpenMS version and output; see [API notes](../references/pyopenms.md).
Done when PSM scores, decoy mapping and thresholded results are reviewed.
