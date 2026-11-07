# Single-amplicon Cas9 or HDR
~~~bash
CRISPResso --fastq_r1 sample_R1.fastq.gz --fastq_r2 sample_R2.fastq.gz --amplicon_seq "$AMPLICON" --guide_seq "$GUIDE" --quantification_window_size 1 --quantification_window_center -3 --min_average_read_quality 30 --output_folder sample_results --name sample_id
~~~
Set AMPLICON to the primer-trimmed reference and GUIDE to the protospacer without PAM. For HDR add the expected edited amplicon with --expected_hdr_amplicon_seq. The Q30 filter changes retained reads and potentially measured editing; report aligned and input reads. Inspect window geometry, allele table and mapping statistics; see [window and failure modes](../references/interpretation.md).
Done when outputs and retained fraction are reviewed with the intended editing window.
