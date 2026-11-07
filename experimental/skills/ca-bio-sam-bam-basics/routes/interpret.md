# Interpret alignment fields

```bash
samtools flags 99
```

Use `references/fields-and-cigar.md` for SAM fields, FLAG bits, CIGAR and coordinate conversion. Secondary alignments represent alternative placements; supplementary alignments represent split pieces needed for SV and fusion analyses. For MAPQ choose `references/mapq-by-aligner.md`; for optional tags and program-chain provenance use `references/tags-and-provenance.md`.

Done when the requested fields or flags are interpreted in the context of the aligner and downstream task.
