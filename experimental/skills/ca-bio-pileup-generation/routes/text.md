# Text pileup

```bash
samtools mpileup -f reference.fa input.bam > pileup.txt
```

Index the FASTA with `samtools faidx reference.fa`; index the BAM for `-r`. Use `-r chr1:1000000-2000000` for an indexed region or `-l targets.bed` for targets. Multiple BAMs may follow. Samtools defaults cap depth at 8000 and apply BAQ when a reference is supplied; use `-d 0` for uncapped depth when appropriate. For assay-specific filters, text encoding, depth, BAQ and overlap caveats read `references/text-and-filters.md` and `references/options-errors.md`.

Done when the requested positions and filters are represented in the pileup.
