# Search and rescore
~~~bash
MZML="sample.mzML" FASTA="targets.fasta" ENGINE=sage OUT=dda_out bash scripts/dda_search.sh
~~~
Use centroided input, a target-only FASTA and an installed Sage or Comet, OpenMS DecoyDatabase and Percolator. Set ENGINE=comet for the other route. The source example uses trypsin, two missed cleavages, specified modifications, instrument-specific tolerances and 8 threads; inspect and adapt before use. It reuses any nonempty existing target_decoy.fasta without checking provenance, and its final q-filtered files omit headers: use a fresh output directory and preserve original Percolator TSVs. Inspect target/decoy labels, q-value columns and logs. See [CLI and MS-GF+ notes](../references/cli_route.md).
Done when database provenance, search settings, decoy counts and filtered results are reviewed.
