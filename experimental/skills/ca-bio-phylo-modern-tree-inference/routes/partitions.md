# Infer with partitions
~~~bash
iqtree3 -s concatenated.fasta -p partitions.nex -m MFP+MERGE -rcluster 10 -B 1000 -bnni -alrt 1000 -T AUTO --seed 12345 --prefix partitioned
~~~
Prepare a partition file matching the actual alignment; -p uses edge-linked proportional branch lengths. The source example contains fixed illustrative gene boundaries that must not be reused blindly. Inspect merged scheme, per-partition information and output.
Done when partition mapping, linkage, selected models and support are reviewed.
