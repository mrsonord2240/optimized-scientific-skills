# Context

Primer-derived sequence does not represent biological sequence and can mask variants. Default ampliconclip matches the 5-prime end; orientation-aware clipping avoids clipping opposite-strand primers. Amplicon reads share coordinates by design: without UMIs skip coordinate deduplication; with UMIs use UMI-aware consensus. --keep-tag retains stale NM/MD; recompute instead. --tolerance defaults to 5; --clipped drops unmatched reads and --fail marks them QCFAIL.
