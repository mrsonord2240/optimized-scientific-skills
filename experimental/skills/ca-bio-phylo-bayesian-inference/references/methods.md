## When to Use Bayesian vs ML

| Factor | ML (modern-tree-inference) | MrBayes | BEAST2 | RevBayes | PhyloBayes-MPI |
|--------|----------------------------|---------|--------|----------|----------------|
| Reach for it when | topology, fast, IC model selection, bootstrap | posterior clade PP on a standard alignment; easiest entry; model averaging | TIME trees, tip-dating, demographics, phylodynamics | the model is not a built-in option; full graphical-model control | DEEP phylogeny with compositional heterogeneity / LBA |
| Output | point tree + bootstrap | posterior + clade PP | dated posterior + HPDs | posterior, custom | posterior under CAT/CAT-GTR |
| Speed | fast | moderate | slow | moderate | brutal (days-weeks, MPI) |
| Multi-run | automatic resampling | `nruns=2` default | run seeds independently | manual | run chains independently |
| Citation | Minh 2020 | Ronquist 2012 | Bouckaert 2019 | Hohna 2016 | Lartillot 2013 |

Default recommendation: start with ML (modern-tree-inference) for the topology, then move to Bayesian when posterior probabilities, model averaging, divergence times, or site-heterogeneous models are specifically needed. Need TIMES -> divergence-dating. Need a coalescent species tree -> species-trees.

## Tool Taxonomy

| Tool | Citation | Mechanism / role | When |
|------|----------|------------------|------|
| MrBayes 3.2 | Ronquist 2012 | NEXUS `lset`/`prset`; native MC3 + 2 runs; live ASDSF | general-purpose first choice for a concatenated alignment; partitioned data; `nst=mixed` model averaging |
| BEAST2 | Bouckaert 2019 | BEAUti -> XML; operator tuning; no MC3 by default | TIME trees, tip-dating, coalescent/birth-death demographics, phylodynamics (-> divergence-dating) |
| RevBayes | Hohna 2016 | Rev probabilistic-graphical-model scripting | a model that is not a checkbox: novel mixtures, custom hierarchical priors, biogeography |
| PhyloBayes-MPI | Lartillot 2013 (CAT: Lartillot 2004) | data-augmentation MCMC, MPI-parallel; CAT/CAT-GTR | deep phylogeny, compositional heterogeneity, LBA; `bpcomp`/`tracecomp` convergence |

Decision tree: standard alignment wanting clade PP, easiest path -> MrBayes (gammadir brlenspr, nruns>=2, gate on ASDSF<0.01). Need divergence TIMES -> BEAST2 (-> divergence-dating). Need a model that is not built-in -> RevBayes. Deep phylogeny with suspected compositional heterogeneity / LBA -> PhyloBayes CAT-GTR (or IQ-TREE C60/PMSF in modern-tree-inference if too big for CAT or fast bootstrap is wanted).

## MCMC Convergence -- the Heart

Almost all Bayesian-phylo failures live here, and the non-obvious core is that scalar-parameter convergence and tree-topology convergence are DIFFERENT problems. Continuous parameters (tree length, gamma shape, rates, frequencies) live in a smooth low-dimensional space that mixes well; topology lives in the vast discrete space of (2n-5)!! unrooted trees where a chain can get stuck on one island, mix perfectly within it, post a beautiful trace and ESS for every scalar, and never visit a competing topological peak. Diagnose BOTH; topology is the harder, more-often-skipped one.

- Question A (scalars converged?): ESS, autocorrelation, PSRF on the `.p`/`.log` traces. Tools: Tracer, coda, MrBayes `sump`.
- Question B (the TREE distribution converged?): between-run split-frequency agreement (ASDSF, max split-freq diff), topological ESS. Tools: MrBayes `sumt` (live ASDSF), PhyloBayes `bpcomp`, RWTY (Warren 2017). Passing A and stopping is the single most common silent error.

What each diagnostic measures:
- ESS (effective sample size) = N / integrated autocorrelation time. Check PER parameter (the likelihood can have ESS=2000 while a nuisance rate has ESS=12). Tracer flags <200 yellow, <100 red.
- PSRF (Gelman-Rubin Rhat) compares between-run to within-run variance across >=2 runs; PSRF -> 1.00 at convergence. Requires multiple runs. Scalar-only (Question A).
- ASDSF (average standard deviation of split frequencies) = THE between-run topology metric in MrBayes, computed automatically when `nruns>=2`: for each bipartition above a minimum frequency, take the SD of its frequency across runs and average. Converged runs give the same frequency in both -> ASDSF -> 0. Printed live during `mcmc`; `stoprule=yes stopval=0.01` auto-halts on it, which is a topology-only criterion (see below). Companion: the MAX split-freq diff (worst single split) catches one unsettled clade that a low average hides.
- Multiple independent runs are mandatory: PSRF and ASDSF are between-run metrics, undefined with one run, and a single run can look perfectly stationary on a local peak. A single long run is never sufficient evidence of convergence.
- MC3 (Metropolis-coupled MCMC) is MrBayes' defense against tree-space trapping: 1 cold chain (sampled) + heated chains whose flattened posterior crosses valleys, with periodic state swaps. Tuning signal is the chain-swap acceptance rate (target ~20-70%); if near 0, chains too far apart -> lower `temp`; near 100%, too similar -> raise `temp`. BEAST2 and PhyloBayes do not use MC3 by default.

Convergence checklist (all must pass before shipping):
1. Run >= 2 independent runs from different starting trees.
2. ESS > 200 for EVERY parameter (posterior, prior, likelihood, and all model params), not just the likelihood.
3. PSRF ~ 1.00 (<= 1.01) for all scalars (MrBayes `sump`).
4. Trace plots stationary and well-mixed ("fuzzy caterpillar"), independent runs overlay (Tracer).
5. Topology converged: ASDSF < 0.01 (MrBayes) / bpcomp maxdiff < 0.1 (PhyloBayes), and check RWTY tree-space MDS that runs overlap (not separate clouds).
6. Burn-in set from the trace / cumulative split-freq plots, not a dogmatic percentage. If any check fails, run longer; do NOT just increase thinning.

## Priors That Bite

Defaults are not neutral. The branch-length prior is the one that changes conclusions. The classic MrBayes default before 3.2.3 put an i.i.d. Exponential(10) prior on each branch with no control over the SUM, so as taxon number grows the implied prior on total tree length grows and the posterior is pulled toward implausibly long trees with degraded mixing -- and the inflation can feed back into TOPOLOGY and node support, so it is not cosmetic (Brown 2010). The fix is the compound (gamma-)Dirichlet prior (Zhang 2012): a diffuse Gamma on the whole tree length, partitioned among branches by a Dirichlet, decoupling "how long is the tree" from "how is length distributed". It yields posterior tree lengths close to the ML estimate and is robust to its hyperparameters; it became the MrBayes default at 3.2.3+. Rule: never trust a default `exp(10)` brlenspr on datasets with many taxa or long branches -- use `unconstrained:gammadir(...)`; inflated tree length vs an ML reference points first at the branch-length prior. Two further notes: a uniform topology prior is NOT uniform on clades; and Bayes-factor model comparison is only valid under PROPER priors (an unbounded "uninformative" prior leaves the marginal likelihood undefined).

The star-tree paradox and short-internode overconfidence: when the true internal branch is near zero (an effective polytomy), PP does NOT settle toward the uninformative 1/3 among the three resolutions as data accumulate -- it behaves erratically and can drive toward HIGH support for an arbitrary resolution (Lewis 2005; the theory in Yang 2007). So a high PP on a clade subtended by a very short internal branch is exactly where Bayesian support is least trustworthy. Treat near-zero internodes as soft and cross-check against bootstrap. The polytomy-allowing reversible-jump prior of Lewis 2005 is not available in MrBayes 3.2.7a (`prset topologypr` offers only uniform/speciestree/constraints/fixed), so on MrBayes report such a node as an effective polytomy rather than trying to fix it with a prior.

## Per-Method Failure Modes

### Non-Convergence Reported as a Result
**Trigger:** A single short chain, a stationary-looking trace, a point estimate shipped with no ESS / PSRF / ASDSF and no second run.
**Mechanism:** The "result" is an arbitrary draw from a chain that never demonstrably sampled the posterior; the clade PPs are an artifact of where the chain happened to sit.
**Symptom:** No diagnostics reported; or low ESS / PSRF > 1.05 / ASDSF > 0.01 quietly ignored.
**Fix:** Always >= 2 runs; gate on ESS > 200 (all params) AND ASDSF < 0.01 / maxdiff < 0.1; inspect in Tracer/RWTY; never ship without the diagnostics.

### Scalar Convergence Mistaken for Topology Convergence
**Trigger:** Every scalar has ESS > 1000 and PSRF = 1.00, so convergence is declared.
**Mechanism:** The chain is trapped on one tree-space island, mixing perfectly within it while never sampling a competing topology; scalars are necessary but not sufficient.
**Symptom:** Perfect Tracer traces but ASDSF stuck high, or RWTY tree-space MDS shows two separate clouds for the two runs.
**Fix:** Always run a TOPOLOGY diagnostic (ASDSF / bpcomp maxdiff / RWTY MDS), not just Tracer.

### Branch-Length Prior Inflation (Prior Domination)
**Trigger:** Default `exp(10)` brlen prior (old MrBayes) on many taxa or long branches, or a too-strong tree/clock prior on weak data.
**Mechanism:** The i.i.d.-per-branch prior places mass on long total tree length; the posterior converges to the WRONG distribution, and inflated branch lengths can distort topology and support.
**Symptom:** Posterior tree length far exceeds the ML estimate; topology shifts when the prior is changed.
**Fix:** Use `unconstrained:gammadir(...)`; do a prior-sensitivity check (sample the prior alone, vary hyperparameters); compare tree length to an ML reference.

### Harmonic-Mean Model Selection
**Trigger:** Selecting a substitution/clock/partition model by HME (or AICM) Bayes factors.
**Mechanism:** The HME has effectively infinite variance, is unstable, and favors over-parameterized models -> the wrong model is chosen.
**Symptom:** Bayes factors that change on rerun, or a consistent preference for the most complex model.
**Fix:** Stepping-stone or path sampling for marginal likelihoods; interpret with Kass-Raftery 2lnBF thresholds; ensure proper priors.

### PP Overconfidence Under Misspecification
**Trigger:** Reporting PP = 1.0 as proof, especially site-heterogeneous data analyzed under a homogeneous model.
**Mechanism:** A Bayesian analysis conditions on the model; a misspecified model has no way to hedge, so the posterior concentrates with false confidence on LBA-driven clades.
**Symptom:** PP = 1.0 where the bootstrap is 55%, or two conflicting topologies each with high PP across analyses.
**Fix:** Report the model; cross-check PP against bootstrap (sharp disagreement = warning); use CAT-GTR / PMSF for deep data; treat PP = 1.0 as "decisive under THIS model", not truth.

### CAT Non-Convergence at Scale
**Trigger:** Stopping a PhyloBayes CAT run too early because wall-clock budget ran out.
**Mechanism:** CAT chains converge slowly (millions of cycles); maxdiff is still ~0.3 and effsize < 50 when the tree is read out.
**Symptom:** bpcomp maxdiff well above 0.1, tracecomp effsize below 300.
**Fix:** Run long, >= 2 chains, gate on maxdiff < 0.1 AND effsize > 300; budget days-weeks of MPI time; consider removing the fastest sites or recoding.

## Quantitative Thresholds

| Quantity | Threshold | Source |
|----------|-----------|--------|
| ESS per parameter | > 200 (floor 100; 500-1000 for reported point estimates) | Rambaut 2018 |
| PSRF / Rhat | ~ 1.00, accept <= 1.01; > 1.05 = not converged | Ronquist 2012 |
| ASDSF (MrBayes between-run topology) | < 0.01 for a confident result; < 0.05 adequate for hard data | Ronquist 2012 |
| Max split-freq diff | small, consistent with ASDSF < 0.01 | Ronquist 2012 |
| bpcomp maxdiff (PhyloBayes topology) | < 0.1 good; < 0.3 acceptable | Lartillot 2013 |
| tracecomp effsize (PhyloBayes scalars) | > 300 good; > 50 acceptable; rel_diff small (<~0.1) | Lartillot 2013 |
| Burn-in | discard 10-25% (data-driven; up to 50% for slow chains) | Ronquist 2012 |
| MC3 chain-swap acceptance | ~20-70% between adjacent chains; tune `temp` if outside | Ronquist 2012 |
| Bayes factor (2 lnBF) | 2-6 positive, 6-10 strong, > 10 very strong/decisive | Kass & Raftery 1995 |
| Independent runs | minimum 2 (4 for difficult data) | Ronquist 2012 |
| Posterior probability | >= 0.95 nominal but overconfident under misspecification; PP = 1.0 != proof | Suzuki 2002 |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Convergence declared on one run | PSRF/ASDSF undefined with a single run | run nruns >= 2; gate on between-run metrics |
| Perfect ESS but wrong PPs | tree-space trapping not caught by scalars | add a topology diagnostic (ASDSF / bpcomp / RWTY MDS) |
| Posterior tree length >> ML | default `exp(10)` branch-length prior | switch to `unconstrained:gammadir(...)` |
| Model choice flips on rerun | harmonic-mean estimator used | use stepping-stone / path sampling |
| Bayes factor is undefined / arbitrary | an improper "uninformative" prior | make all priors proper before SS/PS |
| PP = 1.0 on an LBA clade at depth | site-heterogeneous data under a homogeneous model | use PhyloBayes CAT-GTR or IQ-TREE PMSF |
| High PP on a near-zero internode | star-tree paradox fabricates resolution | treat the internode as soft; report as an effective polytomy |
| Low ESS after thinning more | thinning discards information, does not raise ESS | run more generations, not higher samplefreq |
| Run stopped early with low ESS | `stoprule=yes` halted on ASDSF (topology only) | `stoprule=no`, or extend with `mcmc append=yes ngen=<new total>` until every ESS > 200 |
| `sump` says "Use the harmonic mean for Bayes factor comparisons" | legacy MrBayes 3.2.7a banner printed with the HME table | ignore it; compare models with `ss` stepping-stone marginal likelihoods |
