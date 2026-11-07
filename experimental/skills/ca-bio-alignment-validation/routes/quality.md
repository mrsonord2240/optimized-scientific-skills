# Alignment quality

```bash
bash scripts/validate_alignment.sh sample.bam > qc_report.txt
```

Exit 0 means PASS or WARN, 1 a FAIL band, 2 unreadable or no primary records. The Python alternative is `python scripts/validate_alignment.py sample.bam > qc_report.txt`; `-n` samples the file head and biases mapping rates. Short-read DNA bands do not apply unchanged to stranded RNA, ATAC, bisulfite or long reads. See `references/metrics.md` for thresholds, individual metrics and troubleshooting.

Done when the report and exit status are interpreted for the assay.
