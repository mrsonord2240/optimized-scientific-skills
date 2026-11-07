# Align two sequences

```bash
python scripts/align.py sequences.fasta alignment.json --mode global --type dna
```

Input must have exactly two FASTA records. Choose `--mode local` for a shared region, or `--mode semiglobal` for free ends on both sequences. Use `--type protein` for BLOSUM62 with source EMBOSS-equivalent gap scores -11/-1. DNA defaults to match 2, mismatch -1, open -10, extend -0.5; `--reverse-complement` compares both query strands for DNA. Read `references/scoring.md` for end-gap direction, BLAST versus EMBOSS conventions and alternative scoring. Report PID2 as aligned residue pairs only; other definitions are in `references/percent-identity.md`. For substitutions, coordinate mapping and export read `references/alignment-export.md`.

Done when the requested alignment and score are written and the mode, orientation, scoring and identity denominator are reported.
