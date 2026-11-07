# sgRNA Library Design for BE Screens

Read when tiling editing-window-positioned spacers across a protein region and annotating target vs bystander edits.

## sgRNA Library Design for BE Screens

**Goal:** Tile editing-window-positioned spacers across a protein region of interest to enable variant scanning.

**Approach:** For each amino acid in the target region, find NGG-adjacent spacers where the SNV-of-interest base falls in editing positions 4-8 with minimal bystander C/A in the same window. Annotate each spacer with the target and bystander editable-base positions (nucleotide level; `find_be_spacers` does not translate to amino acid changes -- pass its output to variant-calling/variant-annotation (VEP) for predicted protein consequence).

```bash
python scripts/find_be_spacers.py --cds cds.fa --target-aa 130 --target-base C --editor BE4max --out spacers.tsv
# or: from find_be_spacers import find_be_spacers  (function; returns a DataFrame)
```

`scripts/find_be_spacers.py` finds sgRNAs that place `target_base` in the editor-specific window at `target_aa`. Editors: BE3, BE4max, eA3A-BE3, ABE7.10, ABE8.20, ABE8e, evoCDA-BE (an unknown editor raises `ValueError`). Returns spacer, strand, spacer_start, target_positions, bystander_positions, n_bystanders, sorted by n_bystanders; an empty frame with the same columns if no PAM-adjacent window holds `target_base`. Input is upper-cased before the case-sensitive PAM search; reverse-strand spacers are converted back to forward-CDS coordinates before the codon-overlap test.

**Decision rule:** Select spacers with target_positions != empty AND n_bystanders minimized. For variant-by-variant scanning, accept up to 1-2 bystanders if biology of those positions is interpretable; flag for downstream variant attribution.
