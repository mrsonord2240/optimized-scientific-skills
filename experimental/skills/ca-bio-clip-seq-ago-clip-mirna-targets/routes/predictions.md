# Project prediction coordinates
~~~bash
python scripts/targetscan_sites_to_bed12.py --sites targetscan8.normalized.tsv --utr-map utr-map.tsv --assembly GRCh38 --annotation-release GENCODEv49 --targetscan-release 8.0 --output targetscan8.GRCh38.bed12 --manifest targetscan8.GRCh38.manifest.json
bedtools intersect -split -s -wa -wb -a ago_peaks.GRCh38.bed -b targetscan8.GRCh38.bed12 > ago_targetscan.tsv
~~~
Done when version-matched transcript coordinates have been projected, strand and assembly verified, and overlap reported as AGO-supported prediction. Use new output paths: the converter overwrites existing files. The source does not supply a complete seed-scan or peak-calling pipeline. A standard AGO peak does not identify the miRNA; filter candidates with matched expression and label inferred assignments. Separate prediction-only candidates and direct chimeras; do not discard direct noncanonical evidence. See [coordinate contract](../references/workflows-and-commands.md).
