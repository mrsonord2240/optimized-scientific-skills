## Model Comparison Done Right

Bayesian model comparison compares MARGINAL likelihoods (the data probability integrated over all parameters under the model's priors) via Bayes factors. The marginal likelihood is a hard high-dimensional integral. The correct estimators are stepping-stone (Xie 2011) and path sampling, which sample a series of power posteriors interpolating prior (beta=0) to posterior (beta=1); stepping-stone is more accurate per step and is the default recommendation. Baele 2012 showed by simulation and empirically that both substantially outperform the harmonic-mean estimator.

The harmonic-mean estimator (HME) is discredited and must never select a model. It is dominated by the smallest likelihoods in the posterior sample (the prior-favored tail the posterior rarely visits), giving it effectively infinite variance -- it does not converge as samples are added, is unstable run-to-run, systematically overstates the marginal likelihood, and favors over-parameterized models. Interpret Bayes factors on the 2*ln(BF) scale (Kass & Raftery 1995): 2-6 positive, 6-10 strong, >10 very strong/decisive, where 2 lnBF = 2*(lnML_1 - lnML_2). BFs require proper priors. The reversible-jump alternative `lset nst=mixed` sidesteps explicit BFs by sampling the substitution model itself (the 203 GTR rate-class groupings) and reporting each model's posterior probability -- the principled way to account for substitution-model uncertainty by averaging. Read the output from `sump`: it prints a "Model probabilities above 0.050" table, one row per sampled `gtrsubmodel[...]` rate-class grouping with its posterior probability (also written to `<filename>.mstat`; checked on MrBayes 3.2.7a). Report the SPREAD across that table, not the single top model -- a run with several groupings each near 0.1-0.2 says substitution-model uncertainty is real and averaging over it (rather than fixing one GTR submodel) is doing work; a single grouping near 1.0 says the data pin down the rate-class structure and `nst=mixed` mainly confirms `nst=6`.

### Run MrBayes and Verify Convergence Before Trusting the Tree

**Goal:** Produce a Bayesian phylogeny whose clade posteriors are trustworthy, gated on proven convergence rather than a single stationary-looking run.

**Approach:** Set a multiple-hit-correcting model and the compound-Dirichlet branch-length prior, run two MC3 runs, watch live ASDSF, then summarize only after the topology and scalar diagnostics pass.

```
begin mrbayes;
    set seed=12345 swapseed=67890;                      [ record seeds for reproducibility ]
    lset nst=6 rates=invgamma;                          [ GTR+I+G; use nst=mixed for rjMCMC model averaging (below) ]
    prset brlenspr=unconstrained:gammadir(1,0.1,1,1);   [ compound Dirichlet, NOT exp(10): avoids tree-length inflation ]
    mcmc ngen=10000000 nruns=2 nchains=4 temp=0.1       [ 2 runs x (1 cold + 3 heated MC3 chains) ]
         samplefreq=1000 printfreq=1000 diagnfreq=5000
         stoprule=no;                                   [ run the full ngen; do NOT auto-halt on ASDSF ]
    sump burninfrac=0.25 relburnin=yes;                 [ scalar PSRF + ESS, discarding first 25% ]
    sumt burninfrac=0.25 relburnin=yes;                 [ consensus tree + clade PP + ASDSF ]
end;
```

Do not use `stoprule=yes stopval=0.01` as the convergence gate: it checks ASDSF (topology) only, and on easy topologies ASDSF reaches 0 at the first diagnostic while scalars are barely sampled (checked on MrBayes 3.2.7a: halted at generation 5,000 with minimum ESS < 10; a larger `diagnfreq` still halted with an ESS of 42). If `sump` shows any ESS < 200, extend the same runs rather than restarting: `mcmc append=yes ngen=<new total>;` (ngen is the new TOTAL generation count, not the increment), then rerun `sump`/`sumt`.

After the run: confirm `sump` PSRF ~ 1.00 and ESS > 200 for every parameter, `sumt` ASDSF < 0.01 with a small max split-freq diff, and check tree-space convergence with RWTY (`analyze.rwty(list(run1=..., run2=...), burnin=round(0.25 * length(run1$trees)))`, then `makeplot.treespace`; rwty's `burnin` is a number of trees, not a percentage). Distrust any clade where PP is high but the bootstrap or concordance factor (modern-tree-inference) is low.

Prior-only run for the prior-sensitivity check: same `lset`/`prset`, then `mcmc data=no ngen=... filename=prioronly;` and compare its `sump` (e.g. tree length) with the data run.

### Compare Two Models by Stepping-Stone

**Goal:** Select between two models (e.g. GTR+G vs GTR+I+G, or partition schemes) by marginal likelihood, never by the harmonic mean.

**Approach:** Estimate each model's marginal log-likelihood by stepping-stone sampling, then convert the difference to a Kass-Raftery Bayes factor.

```
[ run once per model on the same proper priors; nsteps = power-posterior stones ]
ss ngen=1000000 nsteps=50 samplefreq=100 diagnfreq=1000;
[ MrBayes prints the stepping-stone marginal log-likelihood; record it per model ]
[ 2 lnBF = 2 * (lnML_model1 - lnML_model2); interpret on the Kass-Raftery scale ]
```

Samples per step = ngen / (nsteps + 1) / samplefreq (MrBayes prints "N steps will be used with G generations (S samples) within each step"); keep S at roughly 30-50 or more, so raise `ngen` when you raise `samplefreq`. Defaults on 3.2.7a: `alpha=0.4` (beta values from quantiles of Beta(alpha,1)), `burninss=-1` (one step's worth of generations discarded before the first step).

BEAST2: install the `MODEL_SELECTION` package and run `PathSampler` (path sampling / stepping-stone); RevBayes: `powerPosterior()` + `steppingStoneSampler()`. All require proper priors, or the marginal likelihood is undefined.

### Minimal BEAST2 Run

**Goal:** Get from a BEAUti-built XML to a sampled posterior on the command line (checked on BEAST v2.7.7).

```bash
beast -seed 12345 -overwrite my_analysis.xml   # writes my_analysis.log and my_analysis.trees
```

`beast` needs an input XML, which BEAUti (bundled with BEAST2) builds interactively from an alignment plus a site model, clock model, and tree prior -- there is no command-line-only path to a first XML. After the run, treat the `.log` in Tracer exactly like a MrBayes `.p` file (ESS, PSRF across independent seeds) and summarize the `.trees` file with TreeAnnotator. BEAST2 does not run multiple chains itself: launch the same XML under two or more `-seed` values and combine with LogCombiner before checking convergence.

RevBayes has no equivalent minimal command here: choosing RevBayes over the other three tools already means the model is not one of their built-ins, so the Rev script is necessarily bespoke to that model rather than a fill-in-the-blank template. Start from the official tutorials at revbayes.github.io for the specific model (e.g. partitioned GTR+G, or a custom hierarchical prior), and apply the same convergence gate as MrBayes (`powerPosterior()`/`steppingStoneSampler()` for marginal likelihoods; run >= 2 independent chains and check ESS/PSRF, since RevBayes has no built-in ASDSF-style topology diagnostic).

## Site-Heterogeneous CAT Models for Deep Phylogeny

Standard models (GTR+G, LG/WAG+G) assume all sites share one set of equilibrium frequencies -- they are site-HOMOGENEOUS. Real proteins are not: a buried hydrophobic site and a surface charged site have different profiles. At DEEP timescales this misleads, because saturated sites convergently acquire similar compositions in unrelated lineages and a homogeneous model misreads convergent composition as shared ancestry -> long-branch attraction with HIGH PP (the overconfidence-under-misspecification mechanism made concrete). The answer is CAT / CAT-GTR (Lartillot 2004): an infinite-mixture (Dirichlet-process) model assigning sites to an unknown number of categories, each with its own amino-acid frequency profile, all inferred from the data; CAT-GTR adds one shared GTR exchangeability matrix (the PhyloBayes default with 4-category gamma). It is far more robust to LBA at depth than homogeneous models, at the cost of slow convergence (millions of cycles, days-weeks), which is why PhyloBayes-MPI parallelizes it (Lartillot 2013).

CAT chains converge slowly, so diagnostics are mandatory (run >= 2 chains with different names): `bpcomp -x <burnin> <every> chain1 chain2` reports maxdiff (the largest bipartition-frequency discrepancy across chains, the topology metric) and writes a pooled consensus; `tracecomp -x <burnin> chain1 chain2` reports per-statistic effsize and rel_diff (a PSRF analog). Thresholds: maxdiff < 0.1 AND effsize > 300 = converged; maxdiff < 0.3 AND effsize > 50 = acceptable; otherwise keep running. The ML treatment of the same artifact is IQ-TREE C60/PMSF (Wang 2018): fit the C60 profile mixture once, fix each site's posterior-mean profile, then do fast tree search and bootstrap. CAT-GTR gives a posterior (slow, MPI); PMSF gives bootstrap (fast, scales to huge matrices) -- different output, same disease. Routes to modern-tree-inference.

```bash
mpirun -np 8 pb_mpi -d alignment.phy -cat -gtr -dgam 4 chain1   # CAT-GTR + discrete gamma
mpirun -np 8 pb_mpi -d alignment.phy -cat -gtr -dgam 4 chain2   # second independent chain
bpcomp -x 1000 10 chain1 chain2     # maxdiff (topology); writes bpcomp.con.tre
tracecomp -x 1000 chain1 chain2     # effsize + rel_diff (scalars)
```
