# Align homologous sequences

```bash
mafft --localpair --maxiterate 1000 sequences.fasta > aligned.fasta
```

For fewer than 200 sequences the source defaults to L-INS-i; for full-length similar sequences use G-INS-i and for conserved blocks separated by unalignable regions E-INS-i. For larger inputs use FFT-NS-2, MUSCLE5 super5 or ClustalOmega. Commands and parameter choices: `references/aligners.md`. Screen possible non-homology and orientation before alignment; `references/validation.md` includes the source screen and limitations. If `--adjustdirection` reverses a sequence, MAFFT prefixes its ID with `_R_`; reconcile IDs before metadata joins.

Done when the alignment is written and sequence count, homology, orientation, gaps and plausible conserved domains have been checked.
