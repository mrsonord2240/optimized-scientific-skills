# CRISPResso example commands

The source example script contains synthetic sequences and hard-coded input filenames. Adapt the flags to the actual sample and confirm the installed CRISPResso2 version before running:

```bash
CRISPResso --fastq_r1 sample_R1.fastq.gz --fastq_r2 sample_R2.fastq.gz --amplicon_seq AMPLICON --guide_seq GUIDE --base_editor_output --conversion_nuc_from C --conversion_nuc_to T -o results -n sample
```

The source example also shows ABE, prime-editing and batch command variations. Use the relevant complete workflow in `bio-crispr-screens-crispresso-editing` or `bio-crispr-screens-prime-editing-screens` for those tasks.
