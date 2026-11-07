# Multi-scale scPrinter scores

```bash
python scripts/scprinter_footprint.py --fragments fragments.tsv.gz --fasta genome.fa --gtf genes.gtf --blacklist blacklist.bed --regions motifs.bed --outdir sc_out --shift 0,0
```

Use a fresh output directory: the script removes an existing `printer.h5ad` there. Fragments must be bgzip- and tabix-indexed. The shift denotes the shift already applied: source default `0,0` for raw fragments and `4,-5` for Cell Ranger fragments. Add `--groups cells.tsv` for barcode-tab-group pseudobulks; default device is `cuda:0`. First use may download model files. See `references/scprinter.md` for fragment creation, modes, per-cluster depth limits and seq2PRINT training caveats.

Done when `footprints.npz` and `center_by_mode.tsv` are produced and bound/unbound controls and depth are evaluated.
