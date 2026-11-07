# Quantify existing amplicon reads
~~~bash
CRISPResso --fastq_r1 pe_sample.fq.gz --amplicon_seq "$AMPLICON" --guide_seq "$SPACER" --prime_editing_pegRNA_spacer_seq "$SPACER" --prime_editing_pegRNA_extension_seq "$EXTENSION" --prime_editing_pegRNA_scaffold_seq "$SCAFFOLD" --quantification_window_size 25 --output_folder pe_results --name sample_id
~~~
Supply existing verified sequences and input reads. The source specifies extension as RTT reverse complement followed by PBS reverse complement; incorrect ordering can silently misclassify outcomes. Inspect reference, prime-edited and scaffold-incorporated rows and alignment; use canonical bio-crispr-screens-crispresso-editing for full output interpretation.
Done when intended, scaffold and indel outcomes and the configured window are reviewed.
