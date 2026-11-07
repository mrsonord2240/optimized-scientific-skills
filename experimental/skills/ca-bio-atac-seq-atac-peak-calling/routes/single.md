# Single sample and fragment-specific peaks

```bash
macs3 callpeak -t sample.bam -f BAM -g EFFECTIVE_GENOME_SIZE -n sample --nomodel --shift -75 --extsize 150 --keep-dup all -q 0.05
```

Filter peaks with the matched blacklist as shown in `references/other-calls.md`. With one biological sample do not claim replicate reproducibility. For paired-fragment BAMPE omit shift and extension; for sub-100 bp fragments, broad domains, blacklist and signal-track commands use `references/other-calls.md`. Avoid double shifting previously Tn5-shifted input; inspect the aligner settings.

Done when the filtered peak file is written and its caller, threshold and lack of replicate IDR are stated.
