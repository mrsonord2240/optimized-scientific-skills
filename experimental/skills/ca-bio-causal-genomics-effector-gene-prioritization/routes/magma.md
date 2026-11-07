# MAGMA gene association

```bash
bash scripts/magma_genebased.sh
```

First edit the script's input paths, reference and output prefix. Default window is 35 kb upstream and 10 kb downstream; compare 0+0 and 50+50 when appropriate. Match LD ancestry and genome build. The wrapper skips gene-set regression below 200 genes; see `references/magma-gene-based.md`.

Done when gene-based results, window and multiplicity threshold are reported.
