# LDAK SumHer (reference for heritability-partitioning SKILL.md)

## LDAK SumHer Pipeline

**Goal:** Alternative h2 and functional enrichment estimate using the LDAK-Thin model for reconciliation with LDSC.

**Approach:** Reformat sumstats to LDAK input -> run `ldak --sum-hers` against the pre-computed LDAK-Thin tagging file -> compare to LDSC.

Run it as `bash examples/ldak_sumher.sh <gwas_ldak.txt> <trait_prefix>`. Step 1 is `ldak --sum-hers` against the pre-computed LDAK-Thin tagging file; steps 2a/2b are the two-step BaselineLD partition (`--calc-tagging` builds the annotated tagging file with `--power -.25`, `--annotation-number 86`, `--annotation-prefix`; then `--sum-hers` against it, no annotation flags). The script's header comments carry the LDAK input header (`Predictor A1 A2 n Z`) and the BaselineLD download and flag notes.

Pre-computed tagging files exist for GBR (HapMap reference); other ancestries require building tagging file via `--calc-tagging`. LDAK SumHer outputs h2 estimate, per-category h2 share, and enrichment with Z-scores.

## Install

```bash
# LDAK 6+
wget https://raw.githubusercontent.com/dougspeed/LDAK/main/ldak6.3.linux
chmod +x ldak6.3.linux
```

## Common Errors

| Error / symptom | Cause | Solution |
|---|---|---|
| LDAK tagging file: build mismatch | Using hg19 tagging on hg38 GWAS sumstats | Tagging files are build-specific; download matched build |
