# Specialist cleaning and splitting

```bash
divvier -divvy input.fasta
```

Divvier splits ambiguous columns; output may have more columns. For BMGE entropy filtering, use the installed-version flags in `references/bmge-divvier.md`: BMGE 1.12 `-h` and 2.0 `-e` differ, and a successful exit can still write no alignment. For per-residue contamination masking, legacy Gblocks or PhyIN second-pass filtering read `references/specialist-trimmers.md`. Do not apply HMMcleaner to fewer than 15 sequences under the source guidance.

Done when the requested cleaned or split output is checked, with method and impact reported.
