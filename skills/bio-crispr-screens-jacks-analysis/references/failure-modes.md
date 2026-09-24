# JACKS failure modes

## Failure Modes

### Efficacy collapsed near zero for all guides

**Trigger:** Screen used a chemistry the model doesn't support (e.g., CRISPRi screen analyzed with JACKS defaults).
**Mechanism:** CRISPRi efficacy is fundamentally different from Cas9-KO efficacy; the Gaussian efficacy prior (mean 1, variance 1) assumes Cas9-KO-like guides.
**Symptom:** Median efficacy <0.2; almost no significant gene effects.
**Fix:** Do not reuse Cas9-KO efficacies or `--reffile` priors. If a CRISPRi/a reference panel with the same library exists, run JACKS on it and pass its grna file via `--reffile`. Neither the CLI nor `runJACKS()` exposes the prior hyperparameters; they are keyword arguments of `jacks.infer.inferJACKSGene` (`mu0_x`, `var0_x`, `mu0_w`, `var0_w`, `tau_prior_strength`), which `inferJACKS` looks up at call time, so they can be overridden for one run:

```python
import functools
import jacks.infer
from jacks.jacks_io import runJACKS

_default_gene_fit = jacks.infer.inferJACKSGene
# Example: a weaker efficacy prior (larger var0_x) for a chemistry whose efficacies are not centred on 1.
# Choose values from a matched reference run, not by tuning until hits appear.
jacks.infer.inferJACKSGene = functools.partial(_default_gene_fit, mu0_x=1.0, var0_x=4.0)
try:
    runJACKS('counts.txt', 'replicatemap.txt', 'guidemap.txt', outprefix='jacks_crispri', ctrl_sample_hdr='Control')
finally:
    jacks.infer.inferJACKSGene = _default_gene_fit   # restore for later runs in this session
```

This relies on JACKS 0.2 internals, not a public API; re-check `help(jacks.infer.inferJACKSGene)` on other versions.

### Cross-cell-line efficacy disagreement

**Trigger:** Pooling screens across cell lines with very different Cas9 expression / chromatin / fitness baselines.
**Mechanism:** Efficacy depends on Cas9 expression and chromatin accessibility; sharing across lines averages real per-line differences.
**Symptom:** Per-line gene effects look noisier than per-line MAGeCK results.
**Fix:** Use Chronos for multi-cell-line screens with screen-quality modeling; reserve JACKS for screens with matched chemistry + cell type / culture conditions. On a multi-cell-line panel, fit per-cell-line gene effects with shared efficacy and pool effects across lines downstream (meta-analysis), not inside JACKS.

### MCMC / variational convergence failure

**Trigger:** A gene's variational updates hit the iteration cap before the lower bound settles (few guides, conflicting guides, very noisy replicates).
**Mechanism:** Each gene is fitted for at most `n_iter=50` updates and stops early once the lower bound changes by < `tol=0.1`; a gene still moving at iteration 50 is returned as-is.
**Symptom:** Not run-to-run variation -- inference is deterministic, so reruns on identical input give identical effects. Instead, effects shift when you change `n_iter`, and the DEBUG log (`Iter 50/50 ...`) shows genes ending at the cap.
**Fix:** The CLI and `runJACKS()` expose no iteration argument. Refit with a higher cap and compare: `jacks.jacks_io.inferJACKS = functools.partial(jacks.infer.inferJACKS, n_iter=500)` before `runJACKS`, then restore it; genes whose effect/std changes sign or crosses -2 were not converged. Adding guides per gene or screens also helps.

### sgRNA-to-gene map mismatch

**Trigger:** Guide map and count matrix use different sgRNA naming conventions (e.g. `BRCA1_1` vs `BRCA1.1`).
**Mechanism:** JACKS builds a sgRNA-to-gene dictionary from the map and keeps only guides found in it; unmatched guides are dropped, not set to NaN.
**Symptom:** Genes missing from the output rows entirely (no NaN values appear).
**Fix:** Standardize naming; sanity check `len(jacks_output) == n_genes_expected`.

### Reference efficacy prior from wrong library

**Trigger:** Using DepMap Brunello efficacy as prior for a screen with a custom TKOv3-style library.
**Mechanism:** Per-sgRNA efficacy is sequence-specific; sgRNAs in one library map to different gene contexts than another.
**Symptom:** Usually an immediate exception, `<sgRNA> has no sgrna reference in <reffile>`, because JACKS 0.2 requires every guide in the map to be in the reference. That check is by ID only: if the IDs happen to match (e.g. both renamed `GENE_1`, `GENE_2`) the run succeeds silently with the wrong efficacies.
**Fix:** Match library exactly, by sequence and not only by ID; if no matched reference exists, run without prior.
