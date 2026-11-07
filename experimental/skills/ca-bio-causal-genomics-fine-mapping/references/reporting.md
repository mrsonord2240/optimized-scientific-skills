# Required Reporting Schema for Fine-Mapping
Every locus reported should carry these columns; missing fields are the most common reviewer complaint.

| Column | Description |
|--------|-------------|
| locus_id | Locus identifier (chr:start-end or sentinel rsID) |
| method | susie_rss / FINEMAP / PAINTOR / SuSiEx / SuSiE-inf |
| L_used | `sum(fit$V > 0)` (effective L; not just the cap passed in; susieR sets pruned effects' prior variance to 0) |
| n_credible_sets | Number of returned credible sets at the chosen coverage |
| cs_size | Variants per credible set |
| cs_purity_min / cs_purity_mean | min and mean `fit$sets$purity[,'min.abs.corr']` per set |
| top_pip_snp | Lead variant in each credible set |
| top_pip | Posterior inclusion probability of top_pip_snp |
| lambda_s | `estimate_s_rss` diagnostic for the locus |
| kriging_outlier_count | Count of SNPs with `|z_obs - z_exp| > 3` |
| ld_panel | 1KG-EUR / UKB-EUR / in-sample / TopMed |
| prior_source | uniform / PolyFun-EUR / PolyFun matched-ancestry baseline-LF / manual coding-variant |
| coverage | 0.95 default; 0.9 or 0.8 documented for low-N |
| n_effective | Sample size passed to susie_rss (Neff for case-control) |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| `IBSS algorithm did not converge` warning | L too small OR LD mismatch | Increase L; run `estimate_s_rss`; check ancestry match |
| Credible set contains all SNPs at locus | LD reference matches discovery poorly; `s_hat` > 0.1 | Switch to in-sample LD or stratify reference ancestry |
| PIPs identical to GWAS p-value rank | Effectively no LD information used; check LD matrix orientation | Verify SNP order in z and R match exactly; check for transposed R |
| Negative eigenvalues in LD matrix | Numerical PSD violation from finite-precision storage | Add small ridge: `R <- R + diag(1e-4, nrow(R))`; or use `Matrix::nearPD` |
| `pip` all ~ 1/p (uniform) | Convergence failure OR all effects pruned | Check `fit$converged`; raise L; check Z scale |
| FINEMAP `Error: SNP names do not match` | .z and .ld SNP order differ | Ensure both are sorted identically; pass matched .snp file |
| Coloc.susie returns NULL | One trait has zero credible sets | Verify both fits succeeded; lower coverage to 0.9 if signal is weak |
| Coloc.susie crashes with `data.table` error (`:= is defined for use in j`) | `z`/`R` passed to `susie_rss` without SNP-ID names, so `coloc.susie`'s internal SNP match is empty | Name `z1`/`z2` and set matching `ld_matrix` dimnames before fitting; see the precondition in `references/coloc-susie.md` |
| SuSiEx output empty, or fails inside its own internal PLINK calls (`Error: No variants remaining after --extract`) | Per-population lists misaligned with reference panels, or `--plink=<path>` omitted (required on SuSiEx 1.1.2, not caught at argument-parse time) | Verify `--sst_file`/`--ref_file`/`--ld_file` share the same population order and `--bp` window; pass `--plink=$(command -v plink)` — see `references/susiex-cross-ancestry.md` |
| FINEMAP `Cannot recognize flag` (`--prob-tol`/`--n-iterations`/`--n-convergence`); `.ld` "Expected N SNPs ... encountered only 1"; or no plain `locus.cred` written | Flags renamed in FINEMAP 1.4.2; `plink --r square` defaults to tab-delimited (FINEMAP needs spaces); output is `<prefix>.cred<k>` per causal-count model, not a plain `.cred` | Use `--prob-conv-sss-tol`/`--n-iter`/`--n-conv-sss`; add `--r square spaces`; read `locus.cred1`, `locus.cred2`, ...; see `references/finemap-cli.md` |
| `dap-g` exits with status 1 | Normal — DAP-G returns 1 even on a fully successful run | Judge success by stdout (`Independent association signal clusters`), never the exit code |
| PolyFun priors do not change PIPs | Passed to `prior_variance` instead of `prior_weights` | Read susieR docs; use `prior_weights=` |
