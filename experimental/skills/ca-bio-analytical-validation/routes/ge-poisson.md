# Genome equivalents and Poisson template presence

```bash
python scripts/ge_and_poisson.py --input-ng 10 30 --vaf 0.001 0.0001 --out ge_poisson_results.tsv
```

Several `--input-ng` and `--vaf` values give every combination, one row each: genome equivalents, expected mutant molecules (lambda = GE x VAF), `p_template_present` (P of at least `--min-molecules` templates, default 1), and `ge_required` / `ng_required` for `--target` probability (default 0.95). `--ge-per-ng` overrides 303.03 and must be stated. The script refuses a VAF above 1 (a percent).

- Lambda of about 3 (2.996) gives 95% presence: about 29,957 GE (98.9 ng) for one 1e-4 variant, about 2,996 GE (9.9 ng) at 1e-3. At 0.1% VAF on 1,000 GE, lambda is 1 and the template is absent 37% of the time whatever the sequencer does.
- The result is sampling only. It excludes extraction recovery, consensus depth, background error and the calling rule; quote it as "template present", not "detected".
- If the request says detection needs more than one supporting molecule, set `--min-molecules`.
- "More depth is not helping": report unique (consensus) coverage and the dedup plateau; past it the ceiling is GE in the tube.

Done when `ge_poisson_results.tsv` exists and the answer states GE, lambda and that it is sampling only.
