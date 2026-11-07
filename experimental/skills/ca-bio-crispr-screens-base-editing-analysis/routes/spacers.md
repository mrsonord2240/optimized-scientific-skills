# Spacer annotation

```bash
python scripts/find_be_spacers.py --cds cds.fa --target-aa 130 --target-base C --editor BE4max --out spacers.tsv
```

Replace target and editor with the actual input. Inspect target_positions: the helper can include rows without the requested target position. It reports nucleotide positions, not protein consequences; see `references/library-design.md` and `references/editor-selection.md`.

Done when target-containing rows and bystander limitations are reviewed.
