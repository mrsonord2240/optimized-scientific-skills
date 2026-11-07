# Selection and codon-aware processing

```bash
t_coffee input.fasta -mode mcoffee -output fasta_aln -outfile aligned.fasta
```

Avoid aggressive column filtering for dN/dS. For the subsequent TCS confidence scoring, thresholded column selection and MACSE marker conversion read `references/selection-analysis-workflow.md`; its example removes columns, so do not use it when the request requires masking without removal. Preserve codon boundaries; MACSE `trimAlignment` is an option. For confidence assessment use `bio-alignment-multiple`.

Done when the chosen confidence treatment and codon frame are verified and downstream input matches the requested masking or filtering.
