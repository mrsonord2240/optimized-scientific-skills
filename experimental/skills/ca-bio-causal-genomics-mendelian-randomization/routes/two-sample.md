# Two-sample MR

```bash
Rscript scripts/twosample_workflow.R --exposure exposure.tsv --outcome outcome.tsv --outdir mr_out --bfile 1kg_EUR/EUR
```

Input TSVs need SNP, BETA, SE, A1, A2, EAF, P, N. Match ancestry and build. The wrapper filters exposure P < 5e-8 and F >= 10, clumps at r2 0.001 within 10 Mb, harmonizes and writes IVW, Egger, median, mode and sensitivity tables. Use --no-clump only for verified independent instruments. Check sample overlap and correct it separately as in `references/overlap-and-bias.md`; MR-RAPS alone does not fix overlap bias. Steiger is a heuristic and may need binary-trait correlations.

Done when estimates, uncertainty, instrument strength, overlap and sensitivity limitations are reported.
