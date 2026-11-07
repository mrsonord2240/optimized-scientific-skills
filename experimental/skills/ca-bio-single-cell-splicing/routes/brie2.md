# BRIE2 per-cell PSI

```bash
brie-count -a splicing_events.gff3 -S sample_list.tsv -o brie_counts/ -p 16
```

Plate input is BAM path and cell ID per row, sorted and indexed. Run quantification and inspect default event filters using `references/brie2.md`; droplet mode uses different inputs and remains subject to the chemistry gate. `scripts/sc_splicing_brie2.py` provides reusable helper functions, not a CLI workflow. Treat cell-level tests as exploratory if donors are ignored.

Done when requested PSI and uncertainty outputs are available and filtering is reported.
