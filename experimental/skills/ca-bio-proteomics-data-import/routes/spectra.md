# Inspect spectra
~~~bash
python scripts/inspect_mzml.py sample.mzML
~~~
Input may be mzML or mzXML; the script selects the matching pyOpenMS reader and prints MS1/MS2 and missing precursor/offset counts. Missing isolation offsets are unknown, not zero-width windows.
Done when format, scan counts and missing metadata are reviewed.
