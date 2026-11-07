# Kinnex single-cell preprocessing

```bash
skera split raw_kinnex.bam kit_adapters.fasta segmented.bam
```

Obtain the exact kit adapter FASTA, then run lima and Iso-Seq refinement as in `references/single-cell.md`. The source only tested a synthetic array; barcode assignment is handled by `bio-single-cell-splicing`.

Done when segmented reads are ready for the next library-specific stages.
