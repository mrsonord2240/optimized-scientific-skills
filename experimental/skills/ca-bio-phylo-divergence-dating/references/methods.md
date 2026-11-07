## Clock Model Selection

The clock governs how substitution rate varies across branches; choosing it wrong biases dates and misstates their uncertainty. Choose by testing clocklikeness, not by defaulting either way.

| Model | Assumption | When | Diagnostic |
|-------|------------|------|------------|
| Strict | one rate for the whole tree | clocklike data: intraspecific, or short-timescale viral, when a clock test does not reject | most efficient; tightest justified CIs |
| UCLN (uncorrelated lognormal) | each branch rate drawn independently from a lognormal | field default for multi-species data with rate variation | the `ucld.stdev` / coefficient-of-variation diagnostic |
| UCED (uncorrelated exponential) | branch rates drawn from an exponential | larger, less-Gaussian rate swings | available; rarely the first choice |
| Autocorrelated (ACLN) | descendant rate centered on parent (Brownian log-rate; Thorne et al. 1998) | deep trees where rate is heritable (generation time, metabolism); MCMCTree default | rate-variance `sigma2` (MCMCTree) |
| Random local clocks | a few inferred, discrete rate-shift points | episodic / clade-specific rate shifts (Drummond and Suchard 2010) | estimates where and how many shifts |

Relaxed-clock work was established by Drummond et al. 2006 (uncorrelated relaxed clocks, "dating with confidence"). The single most useful BEAST2 relaxed-clock diagnostic is the coefficient of variation (CoV) of branch rates, derived from `ucld.stdev`: a CoV posterior abutting 0 (`ucld.stdev` near 0) means rates are effectively constant and a strict clock suffices (gain precision by simplifying); a CoV clearly above 0 with 0 excluded means the relaxed clock is doing necessary work and a strict clock would be falsely precise. If the `ucld.stdev` posterior just recovers its prior, the data cannot indicate how clocklike the lineages are -- report that. Use the CoV for a quick read, but decide strict-vs-relaxed formally by marginal-likelihood comparison (path sampling / stepping-stone; in BEAST2 the MODEL_SELECTION package's PathSampler -- see bayesian-inference for the estimator and Bayes-factor scale).

## Calibration Strategy

Calibrations are the dominant input. The bedrock rule: a fossil is a MINIMUM, not a point -- a clade is at least as old as a fossil assigned to it, and the true divergence is older, so a near-delta prior on a fossil age forces a guaranteed-too-young, falsely precise node.

| Strategy | Encodes | When | Pitfall |
|----------|---------|------|---------|
| Lognormal (offset) node prior | offset = fossil minimum; true age = min + a modest+ gap | one well-justified fossil, modest gap expected | mean/SD chosen by feel sets the answer |
| Exponential (offset) node prior | firm minimum, agnostic about gap size | good minimum, weak idea of the maximum | long tail can pull the node very old |
| Uniform + soft bounds | min from fossil, max from absence/strat, both leaky (Yang and Rannala 2006) | a defensible minimum AND maximum | hard bounds (no tail) over-dictate |
| Total-evidence / tip dating | fossils as dated, morphologically scored tips (Ronquist et al. 2012) | morphology available; want data-driven fossil placement | the morphological clock is shaky |
| Fossilized birth-death (FBD) | all fossils as samples of one diversification process (Heath et al. 2014) | several fossils; want coherent calibration | needs lambda/mu/psi/rho; sampled-ancestor handling |
| Tip dates (sampling times) | calibration from collection dates | measurably-evolving populations: viruses, ancient DNA | requires a verified temporal-signal check first |

Soft bounds (Yang and Rannala 2006) make a bound a quantile, not a wall: a small canonical 0.025 tail of probability is allowed beyond each soft min/max so one bad fossil cannot dominate. Justify every fossil per Parham et al. 2012 (specimen identity, apomorphy-based placement, geochronologic basis, monophyly of the calibration clade, stated reasoning). The total-evidence approach (Ronquist et al. 2012) includes fossils as dated tips scored for morphology so the data, not the user, place each fossil. The fossilized birth-death process (Heath et al. 2014) is the modern coherent tree prior: it models speciation, extinction, and fossil sampling jointly, uses ALL fossils, allows sampled ancestors, and replaces the incoherent practice of multiplying ad hoc node densities. Prefer FBD when several fossils exist; use simple node densities only for one or two transparent constraints. Secondary calibrations (an age borrowed from another study) launder uncertainty -- never use a point, use the full distribution, and flag it.

## Tool Taxonomy

| Tool | Citation | Mechanism / role | When |
|------|----------|------------------|------|
| BEAST2 | Bouckaert et al. 2019 | full hierarchical Bayesian MCMC; FBD, tip-dating, total-evidence, phylodynamic priors | need a posterior, fossils-as-tips, or complex models; 10s-100s taxa |
| MCMCTree / PAML | dos Reis and Yang 2011 | approximate likelihood: a two-step BASEML gradient + Hessian, then MCMC over a Taylor approximation | genome-scale / many loci where full BEAST is infeasible |
| TreePL / r8s | Smith and O'Meara 2012; Sanderson 2002 | penalized-likelihood point estimate; roughness penalty `lambda` set by cross-validation | very large trees (1000s-10000s taxa); accept point estimates + bootstrap CIs |
| LSD2 / treedater | To et al. 2016 | least-squares dating; native tip-dating, very fast | huge tip-dated viral trees; fast rooting and sanity check before a Bayesian run |

MCMCTree's approximate likelihood does the expensive Felsenstein pruning once (computing branch-length MLEs, gradient, and Hessian per partition) rather than every MCMC step, which is what makes thousands of loci tractable; check that BASEML converged because unreliable per-partition branch lengths (saturated/short loci) corrupt the approximation. The infinite-sites plot (dos Reis and Yang 2011) plots posterior CI width against posterior mean age across nodes: under infinite data the relationship becomes linear through the origin, with the residual width set entirely by calibration uncertainty -- if real-data points already hug that line, more sequence data will NOT narrow the dates and the answer is better fossils, not more sites. `scripts/infinite_sites_plot.py` parses `out.txt`'s posterior-mean/HPD table and fits that through-origin line, reporting the slope and R^2 for one run and, given out.txt from a second run with more loci, whether the slope has stabilized:

```
python scripts/infinite_sites_plot.py 1locus/post/out.txt 2loci/post/out.txt
```

TreePL and r8s return one number per node and produce NO uncertainty; a bare PL date is a failure mode, not a result -- bootstrap for CIs and cross-validate `lambda`.

## Temporal Signal Before Tip-Dating

When samples are collected at different times and the population evolves fast enough to accumulate measurable substitutions between dates (a measurably-evolving population: viruses, ancient DNA), the sampling-date differences themselves are the calibrations -- no fossil needed. Verifying temporal signal first is non-negotiable; a short-time-span dataset frequently has no signal yet a Bayesian run will return a tight, entirely prior-driven date.

- **Root-to-tip regression (TempEst, Rambaut et al. 2016):** build a rough ML tree, regress each tip's root-to-tip genetic distance against its sampling date. A genuine signal gives a positive slope (the slope estimates the rate); a negative slope means no usable signal or a wrong root. R^2 is exploratory only (tips are non-independent): treat near-zero / much-below ~0.2 as a red flag, not a formal pass. The x-intercept estimates the TMRCA as a sanity check; large-residual tips are date/contamination/recombination outliers to investigate before dating.
- **Date-randomization test (Duchene et al. 2015):** the formal test. Re-run the dating with tip dates shuffled many times; the real-data clock-rate estimate must fall OUTSIDE the distribution of randomized replicates (no CI overlap). If the real estimate sits inside that cloud, there is no temporal signal and any date is an artifact of the prior.

## Run Prior-Only to Expose the Effective Prior

**Goal:** Determine whether the molecular data inform each calibrated node, or whether the reported posterior is just a truncated prior reflected back.

**Approach:** Run the same model first with no sequence data to obtain the effective (marginal) prior, then with data; compare specified-vs-effective-vs-posterior on every calibrated node.

```bash
# BEAST2: edit the XML so the MCMC samples from the prior only (no likelihood)
# set <run ... sampleFromPrior="true"> (BEAUti: MCMC panel, "Sample From Prior")
beast -seed 1 -prefix prioronly prioronly.xml      # effective prior on every node
beast -seed 1 -prefix withdata  withdata.xml       # full posterior

# MCMCTree (each run in its own directory): usedata=0 gives the effective prior;
# usedata=3 writes out.BV; copy it to in.BV; usedata=2 is the approx-likelihood posterior
mcmctree mcmctree_prior.ctl    # control file has usedata = 0
mcmctree mcmctree_bv.ctl       # control file has usedata = 3  -> out.BV
cp out.BV ../post/in.BV
mcmctree mcmctree_post.ctl     # control file has usedata = 2 in.BV  (PAML 4.10.10: bare "usedata = 2" exits "file name empty.")
```

PAML 4.10 control-file details (checked on 4.10.10): `BDparas` needs a trailing flag (`BDparas = 1 1 0.1 m`; `m` multiplicative, `c` conditional), and the tree file needs an `ntaxa ntree` header line (e.g. `8 1`) before the calibrated Newick.

Multi-locus data (several genes/partitions dated jointly under one tree): set `ndata` to the number of loci (default 1) and concatenate that many PHYLIP blocks into one `seqfile`, each block starting with its own `ntaxa nsites` header line and the same taxon names, one blank line apart -- the old-style `ndata` option documented in PAML's own `scripts/ndata/README.txt` ("the multiple alignments are in one sequence data file, one after another"), with every locus sharing the single treefile. `scripts/mcmctree_setup.py`'s `write_control_file`/`generate_prior_and_posterior_configs` take `ndata` as a keyword (default 1); pass e.g. `ndata=2` for a two-locus concatenated PHYLIP file.

TreeAnnotator writes node ages and HPDs into `[&height_median=...,height_95%_HPD={lo,hi}]` comments; Bio.Phylo leaves them in `clade.comment` (`clade.confidence` is None) and node order differs between trees with different topologies, so match clades by tip set:

```python
import re
from Bio import Phylo

def node_heights(path):
    heights = {}
    for clade in Phylo.read(path, 'nexus').get_nonterminals():
        tips = frozenset(t.name for t in clade.get_terminals())
        med = re.search(r'height_median=([-\d.eE]+)', clade.comment or '')
        hpd = re.search(r'height_95%_HPD=\{([-\d.eE]+),([-\d.eE]+)\}', clade.comment or '')
        heights[tips] = (float(med.group(1)) if med else None, (float(hpd.group(1)), float(hpd.group(2))) if hpd else None)
    return heights

prior = node_heights('prioronly.mcc.tree')   # effective prior summary
post = node_heights('withdata.mcc.tree')     # posterior summary
for tips, (median, hpd) in post.items():
    # if the posterior median and HPD ~ the effective prior, the data did not inform this node
    print(sorted(tips), 'prior', prior.get(tips, 'clade not in prior MCC tree'), 'posterior', median, hpd)
```

For calibrated nodes the prior-only MCC tree may not contain every posterior clade; the per-calibration MRCA ages in the two BEAST `.log` files (`mrca.age(...)`) give the same comparison without depending on MCC topology.

## Check Temporal Signal Before Tip-Dating

**Goal:** Confirm a heterochronous (virus / ancient-DNA) dataset actually contains clock signal before committing to a Bayesian tip-dated run.

**Approach:** Regress root-to-tip distance on sampling date (positive slope, sane intercept, outliers flagged), then run a date-randomization test; only date if the real estimate sits outside the randomized cloud.

```bash
# Build a quick ML tree, then regress root-to-tip distance on sampling date
# (TempEst does the same thing interactively as a GUI, no CLI: beast.community/tempest)
# (bioconda installs IQ-TREE 3 as `iqtree3`/`iqtree`; on IQ-TREE 2.x the binary is `iqtree2`)
iqtree3 -s seqs.fa -m GTR+G -T AUTO --prefix rttree
python scripts/root_to_tip.py rttree.treefile dates.tsv   # dates.tsv: tip <tab> date
# require a POSITIVE slope (the slope estimates the rate); R^2 is exploratory only,
# treat near-zero / << ~0.2 as a red flag, not a formal pass

# Date-randomization test (Duchene et al. 2015), the formal test: shuffle tip dates
# N times, rerun LSD2 on the SAME topology each time, compare the real rate to the cloud
python scripts/date_randomization.py seqs.fa rttree.treefile dates.tsv 20
# real rate OUTSIDE the shuffled range -> temporal signal; INSIDE -> do not report a date

# Fast non-Bayesian tip-dating + CI as a cross-check (LSD2 via IQ-TREE)
iqtree3 -s seqs.fa -m GTR+G --date dates.tsv --date-ci 100 --prefix lsd2   # dates.tsv: tip <tab> date
```

A positive slope is necessary, not sufficient: a dataset sampled over a few months can still give a positive slope (from root placement) with R^2 < 0.1 and a date-randomization cloud that overlaps the real rate, so always run the randomization test.


## Per-Method Failure Modes

### Effective Prior Is Not the Specified Prior, Unchecked
**Trigger:** Several calibration densities plus a tree prior (Yule / birth-death / FBD), run straight to the posterior.
**Mechanism:** Every node must be older than its descendants, so a parent and child density truncate each other, and the tree prior is multiplied in; the marginal prior can look nothing like either typed density (Heled and Drummond 2012; Warnock et al. 2012).
**Symptom:** A tight posterior credible interval is read as "the data nailed it," when it is really the (truncated) prior.
**Fix:** Always run prior-only (`sampleFromPrior="true"` / `usedata=0`); report specified-vs-effective-vs-posterior per node; if posterior ~ effective prior, the data did not inform it.

### Fossil Treated as a Point, Not a Minimum
**Trigger:** A near-delta calibration density centered on a fossil age.
**Mechanism:** A fossil only bounds a clade from below; the true origin is older by an unknown gap, so a point prior forces a guaranteed-too-young age.
**Symptom:** Falsely precise, systematically too-young dates that propagate across the tree.
**Fix:** Use the fossil as the offset/minimum with a backward tail (lognormal/exponential or soft bounds); never a point.

### No Temporal-Signal Check Before Tip-Dating
**Trigger:** Tip-dating a short-time-span virus or ancient-DNA dataset without TempEst + a date-randomization test.
**Mechanism:** With too little accumulated substitution between sampling dates, the data carry no rate information and the prior drives the date.
**Symptom:** Plausible-looking but entirely prior-driven dates; tight HPDs on data that cannot support them.
**Fix:** Root-to-tip regression (positive slope) AND a date-randomization test (real estimate outside the randomized cloud) before any dating run.

### Penalized-Likelihood Point Estimate Reported With No Uncertainty
**Trigger:** A TreePL / r8s date reported as a single number.
**Mechanism:** PL maximizes a penalized likelihood and returns a point; it produces no posterior or CI, and `lambda` controls how clocklike the tree is forced to be.
**Symptom:** "Clade X is 45 Ma" with no interval, and a `lambda` chosen by default rather than cross-validation.
**Fix:** Cross-validate `lambda` (TreePL `prime` + `cv`); bootstrap sites/input trees and re-run to get CIs; never report a bare PL date.

### Over-Tight Calibrations Drive the Posterior
**Trigger:** One or two narrow calibration densities dominating the timescale.
**Mechanism:** A narrow density propagates through the clock and tree prior to set ages everywhere; the data barely move them.
**Symptom:** The posterior barely differs from the prior, and conclusions reverse when a single calibration is tweaked.
**Fix:** Widen / soften bounds, sensitivity-analyze each calibration one at a time, and prefer FBD for coherence across many fossils.

## Quantitative Thresholds

| Quantity | Threshold | Source / rationale |
|----------|-----------|--------------------|
| ESS (posterior, prior, likelihood, every reported parameter) | > 200 | Rambaut et al. 2018 (Tracer); below ~100 unusable |
| Independent MCMC chains | >= 2, posteriors must overlap | convergence cannot be judged from one chain |
| Burn-in discarded | >= 10% (confirm by trace, not rote) | standard practice; verify stationarity |
| Soft-bound tail probability | 0.025 per bound | Yang and Rannala 2006; MCMCTree `pL = pU = 0.025` |
| Root-to-tip R^2 (scripts/root_to_tip.py or TempEst GUI) | exploratory; near-zero / << ~0.2 = weak signal; positive slope mandatory | Rambaut et al. 2016 (tips non-independent, not a formal test) |
| Date-randomization test | real-data rate estimate outside the randomized distribution (no CI overlap) | Duchene et al. 2015 |
| Coefficient of variation of branch rates | abutting 0 -> strict adequate; clearly > 0 (0 excluded) -> relaxed needed | Drummond et al. 2006 |
| Infinite-sites plot (scripts/infinite_sites_plot.py) | points on the linear CI-width-vs-age line -> more sites will not help | dos Reis and Yang 2011 |
| MCMCTree acceptance proportion | ~20-40% (target ~30%); tune `finetune` if outside | PAML practice |
| Smoothing `lambda` (TreePL/r8s) | set by cross-validation, never default | Sanderson 2002; Smith and O'Meara 2012 |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| MCMCTree refuses to run | no root calibration | set `RootAge` in the control file or a calibration on the root node |
| MCMCTree calibration mis-read | `>`/`<` notation on PAML 4.9b-4.9d (parsing bug fixed in 4.9e; 4.10.10 reads `'>0.35<0.55'` as `B(0.35, 0.55)`) | prefer `B()`/`L()`/`U()`, which also expose the tail probabilities |
| `error: BDparas: expect flag ... C for conditional, M for multiplicative` | PAML 4.10 requires a flag after the birth-death parameters | `BDparas = 1 1 0.1 m` |
| `error: file name empty.` with `usedata = 2` | PAML 4.10 wants the BV file name | `usedata = 2 in.BV` (copy `out.BV` from the `usedata = 3` run) |
| `maintree file can have only one tree` | tree file has no `ntaxa ntree` header | first line e.g. `8 1`, then the calibrated Newick |
| Posterior ~ prior for a node age | data uninformative for that node | report it honestly; do not claim the data estimated the age |
| Wide CIs on both rate and root age | rate-time confounding from too few calibrations | add a well-justified calibration; check the rate-vs-root-age correlation |
| Times off by 100x in MCMCTree | unit confusion (no fixed time unit; the user picks one, 100 Myr conventional so ages are O(1)) | keep calibrations and `rgene_gamma` in that same unit; then 0.6 = 60 Ma |
| Dates conflict wildly with independent evidence | unjustified fossil placement | apply the Parham et al. 2012 checklist; recheck the assigned clade |
| Deep dates biased | substitution saturation at fast sites | use slower markers, amino acids, or codon models; remove saturated partitions |
