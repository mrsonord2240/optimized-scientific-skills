# Protocol-declared targeted UMI
~~~bash
python scripts/extract_targeted_umi.py --read1 R1.fastq.gz --read2 R2.fastq.gz --umi-length "$PROTOCOL_DECLARED_UMI_LENGTH" --library-id "$LIBRARY_ID" --protocol-source "$VERSIONED_PROTOCOL_OR_RECORD" --output-fastq targeted.umi.r1.fastq --manifest targeted.umi.manifest.json
~~~
Done when the new FASTQ and manifest bind the library, protocol and declared read-2 UMI length. Pinned upstream prose and script disagree on length; neither is a default. Targeted and total library layouts differ; preserve separate primer handling and confirm the protocol. See [workflow contracts](../references/workflows-and-commands.md).
