# FASTQ files to a guide count table

```bash
mageck count --list-seq library.csv --sample-label Plasmid,Day0,Veh_r1,Veh_r2,Drug_r1,Drug_r2 \
  --fastq Plasmid.fq.gz Day0.fq.gz Veh_r1.fq.gz Veh_r2.fq.gz Drug_r1.fq.gz Drug_r2.fq.gz \
  --norm-method median --output-prefix experiment --trim-5 5
```

`library.csv` columns are read by position, not by name: guide id, sequence, gene (`sgRNA,Sequence,Gene`). Gene in the second column gives 0% mapped. Labels and FASTQs must be in the same order. `--trim-5` takes a base count or `AUTO`, not an adapter sequence; 5 trims the CACCG scaffold.

- Sequence the plasmid pool and keep it as a sample: it is the baseline for the cloning bottleneck.
- Mapping rate below 65-70%: check the trim and the library format before going on.
- `experiment.count.txt` is the raw table every later route reads; `experiment.countsummary.txt` has per-sample Gini and zero counts.
- Cas12a libraries (Inzolia, in4mer) and 10X direct capture: `routes/branches.md`.

Done when `experiment.count.txt` exists. Next: `routes/qc.md`.
