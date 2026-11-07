# Transcript TPM PSI

```python
from scripts.quantify_splicing import run_suppa2_quantification
run_suppa2_quantification('annotation.gtf', 'transcript_tpm.tsv', 'events')
```

The TPM header contains sample names only; `scripts/quantify_splicing.py` provides write_suppa_tpm and filter_reliable_events. SUPPA2 is annotation-bound and version 2.4 requires statsmodels <0.15. See `references/suppa.md`; low-replicate differential testing requires special care in `references/failure-modes-and-errors.md`.

Done when the requested PSI files are generated and missing events are noted.
