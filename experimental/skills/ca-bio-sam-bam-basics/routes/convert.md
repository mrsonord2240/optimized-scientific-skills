# Convert alignment formats

```bash
bash scripts/convert_formats.sh input.bam output.cram reference.fa
```

The output extension selects SAM, BAM or CRAM; the optional third argument supplies the reference for CRAM input or output. Input and output must differ. A SAM lacking @SQ lines needs a reference FAI with samtools `-t`. BAM to CRAM and back can change CIGAR representation and tags; retain BAM when exact representation matters. See `references/cram-reference.md` for full decode checks and embedded-reference options.

Done when output is readable with its matching reference and conversion is validated.
