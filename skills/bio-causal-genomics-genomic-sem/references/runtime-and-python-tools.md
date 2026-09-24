# Computational footprint and Python tool install (reference for genomic-sem SKILL.md)

Runtimes and cluster settings for the genome-wide steps, and the LDSC / MTAG Python installs. Read it when `SKILL.md` "Reference Files" points here. The section names below refer to `SKILL.md`.

## Computational Footprint

| Step | Runtime | Hardware |
|------|---------|----------|
| `ldsc()` multi-trait sampling covariance | minutes | laptop |
| `commonfactor()` / `usermodel()` (no SNP loop) | seconds | laptop |
| `commonfactorGWAS()` over 6-8M SNPs | 4-24h depending on cores | cluster recommended |
| `userGWAS()` over 6-8M SNPs with complex path model | 8-48h | cluster recommended |
| Stratified GenomicSEM with 50+ annotations | 1-3 days | cluster |
| MTAG over 6-8M SNPs | 1-2h | laptop or cluster |

Cluster runs of `commonfactorGWAS()` / `userGWAS()` should use `MPI=TRUE` when submitting via mpirun; GenomicSEM detects the OS internally (via `Sys.info()[['sysname']]`) and selects FORK (Linux/Mac) vs PSOCK (Windows) cluster types automatically -- there is no `Operating=` user argument. On a Mac/Windows workstation, reduce `cores` to the physical-core count to avoid PSOCK fork failures.

## Tool Installation (Python tools)

For Python tools:

```bash
# LDSC (munge_sumstats.py for GenomicSEM input): CBIIT/ldsc, checked on commit 1f09cf0,
# Python 3.9. abdenlab/ldsc-python3 v2.0.0 and belowlab/ldsc v3.0.1 crash on --h2/--rg.
git clone https://github.com/CBIIT/ldsc.git && cd ldsc
micromamba create -n ldsc -c conda-forge -c bioconda python=3.9 bitarray=2 pybedtools=0.10.0 -y
micromamba run -n ldsc pip install numpy==1.21.5 pandas==1.3.3 scipy==1.7.3

# MTAG
git clone https://github.com/JonJala/mtag.git
cd mtag && pip install -r requirements.txt
```
