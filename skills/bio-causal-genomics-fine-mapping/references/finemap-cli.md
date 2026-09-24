# FINEMAP CLI

## FINEMAP CLI Pattern

**Goal:** Independent confirmation via shotgun stochastic search.

**Approach:** Build .z, .ld, and master files; run FINEMAP with `--sss` and parse the .snp and .cred outputs.

Full pipeline (locus extraction, `plink --r square spaces` LD, `.z` and master-file construction, run, report): `examples/finemap_pipeline.sh` (`.z` columns: rsid chromosome position allele1 allele2 maf beta se; `.ld` is a square, space-separated, headerless signed-r matrix — PLINK's `--r square` writes tab-delimited by default, so the `spaces` modifier is required, not optional; master header `z;ld;snp;config;cred;log;n_samples`). Tighter search than the example's defaults:

```bash
finemap --sss --in-files locus.master --n-causal-snps 5 --prob-conv-sss-tol 0.001 --n-iter 100000 --n-conv-sss 1000
```

**Checked on FINEMAP 1.4.2** (`finemap --help`): the flags are `--prob-conv-sss-tol`, `--n-iter`, and `--n-conv-sss` — `--prob-tol`, `--n-iterations`, and `--n-convergence` are not recognized and abort with `Cannot recognize flag`.

Outputs: `locus.snp` (per-variant PIP `prob`, `log10bf`), `locus.config` (top configurations), and one `locus.cred<k>` file per evaluated causal-count model (e.g. `locus.cred1` when the posterior favors exactly one causal SNP) — FINEMAP does not write a plain `locus.cred`; report the highest-posterior `k`'s file.

FINEMAP and SuSiE agree when sparsity holds; disagreement often reveals non-sparse loci that need SuSiE-inf.
