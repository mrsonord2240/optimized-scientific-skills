# Design and annotation review

```bash
gffread genes.gtf --bed | cut -f1-12 > genes.bed12
```

For planned sequencing, evaluate library, read length, depth and replicates using `references/design.md`; intron retention calls require suitable rRNA-depleted libraries. Annotation choice is covered in `references/library.md`.

Done when design limitations and suitable next steps are identified.
