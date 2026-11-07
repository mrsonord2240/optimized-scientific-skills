# Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| Every AFR sample flagged as high het/hom | One global cutoff across ancestries | Infer ancestry, flag outliers within each group |
| Good variants filtered by HWE | Two-sided HWE on a mixed cohort or on cases | Excess het only, within ancestry, controls only |
| Spurious HWE or missingness | Cohort metrics computed before genotype filtering | Set `./.` on low GQ/DP first, then recompute |
| `bcftools gtcheck -G 1` errors | `-G` removed after the 1.10 rewrite | Drop `-G`; without `-g` it cross-checks all samples in the file |
| Low Ti/Tv only in the novel set | False positives concentrate in novel variants | Raise stringency; recheck against dbSNP overlap |
| `plot-vcfstats not found` or no plots | matplotlib missing or not on PATH | `pip install matplotlib`; `which plot-vcfstats` |
| `plot-vcfstats` exits 2 at the PDF step | No pdflatex or tectonic | Install a LaTeX engine for `summary.pdf`; the PNGs are already written |
| peddy `IndexError` in `par_het`, or nonsense | The ~25000-site panel and chrX assume a human genome-wide callset | Use somalier or `gtcheck` on other organisms or targeted panels |
| somalier finds 0 or few sites | `--sites` panel built for another genome build, or wrong for a custom panel | Use the sites file for your build, or `somalier find-sites` |

Done when the failing command has run to completion.
