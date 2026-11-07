## SUPPA2 Workflow

**Goal:** Compute event PSI from transcript TPM without alignment; useful when Salmon/kallisto TPMs already exist.

**Approach:** Generate IOE event definitions from GTF, then aggregate TPMs of transcripts including/excluding each event.

```bash
suppa.py generateEvents -i annotation.gtf -o events -f ioe -e SE SS MX RI FL

# -e takes {SE,SS,MX,RI,FL}; SS emits A5+A3 files, FL emits AF+AL files
for ev in SE A5 A3 MX RI AF AL; do
    suppa.py psiPerEvent -i events_${ev}_strict.ioe -e transcript_tpm.tsv -o psi_${ev}
done
```

Install SUPPA2 with `statsmodels<0.15` (0.15 breaks `suppa.py`). The TPM table is tab-separated with a header of sample names only, no cell above the transcript IDs (`write_suppa_tpm` in `scripts/quantify_splicing.py` writes it; a pandas default header gives "6 expected, 5 given" and no output).

SUPPA2 is annotation-bound: events absent from the GTF cannot be quantified. Whether an event is detected depends entirely on which transcripts the upstream Salmon/kallisto index contains. Use GENCODE comprehensive over basic when SUPPA2 detection sensitivity matters.
