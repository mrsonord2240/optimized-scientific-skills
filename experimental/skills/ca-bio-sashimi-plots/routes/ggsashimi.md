# Grouped overlay

```bash
ggsashimi.py -b groups.tsv -c chr17:43094000-43125000 -o figure -g annotation.gtf --alpha 0.25 --height 3 --width 10 --shrink --fix-y-scale -O 3 -C 3 -P palette.txt -A mean_j -F pdf
```

The headerless TSV has sample, BAM path, and group columns; palette has one color per group. Check all paths, region reads and nonempty PDF: the tool may exit 0 on failure. For validated Python plotting and batch rMATS selection use `scripts/plot_sashimi.py` as documented in `references/ggsashimi.md`. Default -M 1 is inclusive per sample before averaging; higher values can bias labels. mean_j labels are rounded raw junction counts, not depth-normalized. For stranded data choose the correct -s mode; paired and single-end modes differ. If no junction survives, drop --shrink.

Done when the figure exists, arcs and exon tracks align, and group colors and scales are clear.
