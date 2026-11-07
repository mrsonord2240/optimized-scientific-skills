# Summarize label-free data
~~~bash
Rscript scripts/msstats_summarize.R evidence.txt proteinGroups.txt annotation.csv protein_abundance.csv
~~~
The source uses MSstats equalizeMedians, Tukey median polish and MBimpute FALSE; optional fifth argument TRUE enables the AFT censored route. Verify annotations and input row counts. For normalized MaxLFQ ratio extraction:
~~~bash
Rscript scripts/maxlfq_iq.R peptide_long.csv protein_maxlfq.csv
~~~
The CSV needs protein, ion, run and raw intensity. iq first normalizes per run and reports disconnected protein sample sets in the sidecar; disconnected components do not share a scale. Median centering of a protein matrix is a normalizer, not MaxLFQ:
~~~bash
python scripts/normalize_matrix.py median raw_intensities.csv centered.csv
~~~
Inspect missing values and the assumption that most proteins are unchanged.
Done when inputs, summarizer, normalization and disconnected groups are documented.
