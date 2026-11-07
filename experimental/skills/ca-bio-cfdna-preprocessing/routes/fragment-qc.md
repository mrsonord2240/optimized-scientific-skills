# Fragment-size QC

```bash
python scripts/cfdna_cli.py qc aligned.bam qc.json --max-size 600
```

The input is a BAM; the output is a JSON summary. The source counts proper-pair primary PF nonduplicate records with positive template length up to the maximum; reason counts overlap. Empty input yields zero summaries, so inspect `n`. Interpret mode, median, 90–150 bp fraction and fraction over 250 bp in light of chemistry, collection and any size selection. Duplicate flags can themselves bias a quantitative cfDNA sample; coordinate-only deduplication can collapse independent molecules. Avoid blanket MAPQ or length filtering for fragmentomics without considering its bias.

**Done when:** the summary and filtering, chemistry and pre-analytic context are retained.
