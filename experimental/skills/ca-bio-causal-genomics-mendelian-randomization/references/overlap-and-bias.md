# One-Sample vs Two-Sample Bias Direction
| Design | Weak-IV bias direction | Reason |
|--------|------------------------|--------|
| One-sample, F<10 | Toward confounded observational estimate (overestimates causal effect if confounding is in same direction) | Sample correlation between IV-X and IV-Y residuals |
| Two-sample non-overlapping, F<10 | Toward null | Independent samples decouple residuals (Burgess 2011 IJE 40:755) |
| Two-sample with partial overlap | Intermediate; proportional to overlap fraction and z-score correlation | Burgess 2016 Genet Epidemiol 40:597; correction available |

**Operational rule:** Whenever both GWAS came from UK Biobank (or any single biobank), treat the analysis as one-sample-equivalent. MR-RAPS corrects weak-instrument bias in this regime; it does NOT correct the separate sample-overlap/confounding bias that one-sample-equivalent designs also carry, and the two can coexist even when F is well above the one-sample floor (F >= 20) -- a real UKB-on-UKB run at mean F=44.6 with a planted null effect showed naive IVW biased (b=0.18) and MR-RAPS no better (b=0.196), because the bias source was sample-overlap confounding, not weak instruments. Use MR-RAPS when F is borderline; use MRlap (below) or Burgess 2016 overlap-corrected IVW whenever overlap/confounding is suspected, regardless of F. Treating it as "two-sample because separate GWAS files" is a common error and produces overestimates.

MRlap (joint overlap + winner's curse + weak-IV correction): code, decision rule and the h^2 < 0.05 fallback are in `references/mrlap-overlap-correction.md`. Prefer it whenever any sample overlap is suspected, and always for UKB-on-UKB.
