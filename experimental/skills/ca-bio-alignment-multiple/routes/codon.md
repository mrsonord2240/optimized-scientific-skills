# Codon-aware alignment

```bash
mafft --localpair --maxiterate 1000 proteins.fasta > proteins_aligned.fasta
pal2nal.pl proteins_aligned.fasta codons.fasta -output fasta > codons_aligned.fasta
```

Match protein and CDS identifiers and verify translations and reading frames. For a non-standard genetic code pass PAL2NAL `-codontable N`; mixed codes require separate translation by lineage. For frameshifts, PRANK, MACSE, OMM_MACSE or HyPhy workflows use `references/codon-methods.md` and `references/specialist-tools.md`.

Done when the codon alignment is written and IDs, translation and frame are validated.
