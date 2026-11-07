# Convert alignments

```bash
python scripts/alignment_io.py convert input.aln clustal output.fasta fasta
```

Select formats for the actual files. For NEXUS, add `--molecule-type DNA`, `RNA`, or `protein` based on the data. Strict PHYLIP and Clustal can truncate IDs; verify uniqueness after writing. Stockholm conversions may discard annotations. See `references/formats.md` for formats, annotation preservation and downstream tool requirements. Use `inspect` instead of `convert` to obtain a JSON summary of all alignments.

For programmatic alignment construction, directory batches, or Bio.Align objects, use `references/advanced.md`.

To inspect aligned sequences themselves, write FASTA using the conversion command and read the resulting file. For a single sequence or column range, use `routes/slice.md`.

Done when the output is written and its format and identifiers checked.
