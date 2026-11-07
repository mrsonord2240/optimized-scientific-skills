# Library characteristics

```bash
infer_experiment.py -i sample.bam -r genes.bed12 -s 200000 > strand.txt
```

Interpret the dominant strand and map it to the downstream tool's convention. For Picard RNA metrics, gene-body bias, rRNA fractions and annotation choice see `references/library.md`; check low-depth bias with gene-body coverage.

Done when the requested strand, coverage bias or contamination is interpreted.
