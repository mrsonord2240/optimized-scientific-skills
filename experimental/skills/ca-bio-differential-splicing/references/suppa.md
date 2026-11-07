## SUPPA2 Differential Analysis

**Goal:** Quick differential splicing from existing transcript quantifications, useful as a sanity check or pilot.

**Approach:** Generate per-condition PSI files from Salmon TPM, then run `diffSplice` with empirical or classical p-values.

```bash
suppa.py generateEvents -i annotation.gtf -o events -f ioe -e SE SS MX RI

for ev in SE A5 A3 MX RI; do
    suppa.py psiPerEvent -i events_${ev}_strict.ioe -e ctrl_tpm.tsv -o ctrl_${ev}
    suppa.py psiPerEvent -i events_${ev}_strict.ioe -e trt_tpm.tsv -o trt_${ev}

    suppa.py diffSplice \
        -m empirical \
        -gc \
        -i events_${ev}_strict.ioe \
        -p ctrl_${ev}.psi trt_${ev}.psi \
        -e ctrl_tpm.tsv trt_tpm.tsv \
        -o diff_${ev}
done
```

The TPM files need a header of sample names only (no `Name`/`transcript` label); generateEvents on an SE-only annotation writes empty A5/A3/MX/RI files. SUPPA2 2.4 imports `multipletests` from a statsmodels module that statsmodels 0.15 removed: use statsmodels <=0.14 (0.14.6 worked).

**Minimum replicates: n>=4 per group, in either mode.** Measured on a simulation (24 strong + 6 weak planted ΔPSI events among 120 genes, TPM from the same molecules; calls at p<0.05 and |ΔPSI|>0.10):

| n per group | empirical: strong found, false calls / all calls, calls on a null comparison | classical: strong found, smallest p |
|-------------|------|------|
| 2 | 19/24, 6/26, 8 | 0/24, p 0.22 |
| 3 | 21/24, 4/28, 1 | 0/24, p 0.077 |
| 4 | 22/24, 1/24, 0 | 21/24, p 0.027 |
| 5 | 21/24, 0/22, 0 | 24/24, p 0.008 |

- **Classical (Wilcoxon) is not a fallback for n<=3.** An unpaired rank-sum test cannot reach p<0.05 there (BH cannot lower it), so it found nothing at 2v2 and 3v3 (real chrX 2v2: smallest p 0.22).
- **Empirical at n<=3** finds most strong events but 14-23% of its calls are false and it calls events on a null comparison (real chrX 2v2, permuted labels: 24 events; true labels: 36). Use it only to screen candidates and confirm with rMATS/leafcutter/Shiba; for n<=3 prefer those tools.
- ΔPSI = group 2 - group 1.
