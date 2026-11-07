# Direct chimeric assignments
~~~bash
bash scripts/run_chimeric_eclip.sh --reads R1.trim.fq.gz --hyb-db mirna_transcripts_GRCh38_GENCODEv49 --run-id sample --output-dir results
~~~
Done when the named prebuilt Hyb database and reference versions are recorded and consensus sites, targets, support, exclusions and manifest are reviewed. The wrapper requires the pinned Hyb interface and at least two replicates; it excludes unstable or missing assignments. Optional matched expression requires both --expressed-mirnas and --expression-threshold. Counts reflect recovery, not affinity. Verify exact installation and read layout before use; no database is built by the wrapper. See [workflow contracts](../references/workflows-and-commands.md) and [failure modes](../references/failure-modes.md).
