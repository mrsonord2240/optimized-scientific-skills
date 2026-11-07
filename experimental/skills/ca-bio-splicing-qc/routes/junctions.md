# Junction quality report

```bash
python scripts/splicing_qc.py report sample.bam genes.bed12 sample_qc > report.txt
```

The helper checks BED12, shared contigs and spliced reads, skips R plotting unless requested, and reports read-weighted known fraction, saturation and fragment-level support. For high novel fraction use `python scripts/splicing_qc.py novel-qc sample.bam sample_qc.junction.xls --min-overhang 8`. A high novel fraction can reflect biology or artifacts; see `references/junctions.md` and `references/troubleshooting.md`.

Done when saturation, known/novel fractions and junction support are interpreted.
