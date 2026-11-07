# Synthetic replicas of two variantbench problems (2026-09-19)

Everything here is **fabricated**. The originals' data live behind Latch and are not obtainable, so these re-create the *structure* of two problems from `github.com/latchbio/variantbench` (Apache 2.0). Results are **not** benchmarks.bio scores and must never be reported as such.

| dir | replicates | tests |
| --- | --- | --- |
| `fine_mapping_effector/` | `fine_mapping_effector_gene_prioritization` | integrate coding, cis-eQTL and trans evidence; don't follow a single heuristic |
| `cis_mr_bcaa/` | `cis_mr_bcaa` | cis-MR at proxy loci when the target gene has no instruments; avoid pleiotropic genome-wide MR |

Design source was each original's public `notes` field (intent, traps), not our Skills. Ground truth is known by construction.

## What the agent may see

Only the `task` string from `eval.json` and the contents of `data/`. **Never** expose: `graders`, `truth.json`, `generate.py`, `solvers.py`, this README. `graders` hold the answer key; `truth.json` holds the true causal effects.

## Validation

`python solvers.py` in each dir runs a reference solver that should pass and several single-heuristic solvers that must fail, through LatchBio's own graders (`common/grading.py` loads them from `../latch-eval-tools` without the Latch SDK).

- **fine-mapping:** integrated solver passes; traps A, C, D, E, F and "count synonymous/UTR as coding" all fail.
- **cis-MR:** correct solver passes; genome-wide IVW (est 0.236, "YES"), direct-at-BCKDK (no instruments) and wrong-proxy-set all fail.

**This proves the answer key and graders are wired correctly. It does not show the tasks are hard for an LLM.** Only a no-skill pilot run shows that.

## Known limits

- **cis-MR blind spot:** true effects (|θ| < 0.1) are smaller than the ±0.15 estimate tolerance inherited from the original, so an allele-harmonisation sign error is *not* detected by the estimate grader. Verified with a solver that skips harmonisation: it passes. Loci set, pleiotropy trap and YES/NO are the discriminating parts.
- **fine-mapping:** the original leans on recalling real biology (a named gene's role in interferon signalling). Here that is supplied as a `gene_description` / `annotation` column, which makes it easier than the original. Gene symbols are invented (`SYN*`) so nothing can be recalled.
- **cis-MR:** pathway gene names (BCKDK, DBT, PPM1K, BCAT2 ...) are real because the prompt names them; coordinates, genotypes and every statistic are simulated. LD structure is a latent-Gaussian threshold model, not a coalescent simulation.
- **Seed 1** was chosen from seeds 1-6 for cis-MR so the realised estimates sit near the designed truth with no locus individually significant. That is a constraint on the fabricated world, recorded in `generate.py`.
- Estimates are asked as log odds ratios (original left the scale unstated).
- The PLINK reference panel needs `plink`/`plink2` or a bed reader in the agent container. Verified against the independent `bed-reader` package.

## Reproduce

```
cd fine_mapping_effector && python generate.py && python solvers.py
cd ../cis_mr_bcaa && python generate.py && python solvers.py --write-eval && python solvers.py
```
Needs numpy, scipy, pandas (+ `bed-reader` only for the cross-check).
