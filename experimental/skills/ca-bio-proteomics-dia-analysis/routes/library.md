# Search using an existing library
~~~bash
diann --f sample1.mzML --f sample2.mzML --lib spectral_library.parquet --fasta database.fasta --out diann_out/report.parquet --qvalue 0.01 --matrices --mass-acc 15 --mass-acc-ms1 15 --reanalyse --smart-profiling --threads 8
~~~
Supply a compatible library, FASTA and data; replace source example mass accuracies with validated instrument settings. This is a library search with reanalysis and no FASTA-search/predictor flag. Version-specific flags, modification consistency and library coverage must be checked.
Done when the report, log, library provenance and search settings are reviewed.
