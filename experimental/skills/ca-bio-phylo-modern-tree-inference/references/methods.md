## Tool Taxonomy

| Tool | Citation | Role | When |
|------|----------|------|------|
| IQ-TREE2 | Minh 2020 | ML search + ModelFinder + UFBoot2 + SH-aLRT + gCF/sCF + AU test + C60/PMSF, all built in | the default for almost all work |
| RAxML-NG | Kozlov 2019 | ML search, transfer bootstrap, terrace-aware, MPI/checkpointing | very large trees, TBE, precise branch lengths, long HPC runs |
| PhyML | Guindon 2010 | ML search, origin of SH-aLRT | legacy/teaching; SH-aLRT is now in IQ-TREE2 |
| FastTree | Price 2010 | approximate ML, single-pass NNI/SPR | a fast first-pass tree on thousands of sequences; not for final support |

IQ-TREE2 vs RAxML-NG (both hill-climb the same likelihood; they differ in built-in features and scaling, not correctness):

| Need | Use |
|------|-----|
| Model selection, UFBoot2, SH-aLRT, concordance factors, AU test, mixture/PMSF | IQ-TREE2 (built in; nothing else bundles all of this) |
| Very large trees (thousands of taxa), low memory, MPI | RAxML-NG |
| Transfer bootstrap (TBE) for rogue-taxon mega-trees | RAxML-NG (`--bs-metric tbe`) |
| Most precise branch lengths for downstream dating | RAxML-NG |
| Robust checkpoint/restart on long runs | RAxML-NG |

Common production pattern: model-select, compute CFs, and topology-test in IQ-TREE2; do heavy bootstrap on a mega-tree in RAxML-NG with `--bs-metric fbp,tbe`.

## Model Selection

ModelFinder (Kalyaanamoorthy 2017) scores the substitution matrix and the rate-heterogeneity model jointly, including FreeRate categories the old jModelTest/ProtTest generation never tested, and ranks by BIC (the default; k*ln(n) penalty favors simpler models that generalize on large n).

- `-m MFP` -- ModelFinder Plus: test all models by BIC, then search with the winner. The standard. (`-m MF` selects only; avoid `-m TEST`/`-m TESTONLY`, the legacy jModelTest-style limited set.)
- **+G vs +R.** `+G` (discrete Gamma, one alpha) is the workhorse for single short genes. `+R` (FreeRate) freely estimates each category's rate AND weight, capturing the non-Gamma, often multimodal rate distributions of concatenated phylogenomic data; expect `+R3` to `+R6` to win on BIC at scale.
- **The +I+G trap.** The invariant-sites proportion (+I) and the Gamma shape (+G) describe overlapping rate distributions, so the likelihood surface has a flat ridge: the two estimates are individually near-meaningless and start-dependent. Prefer `+R` (its slowest category absorbs near-invariant sites); use `+G` alone unless BIC genuinely demands +I+G.
- **Partition model selection.** `-m MFP+MERGE` fits per-partition models AND greedily merges partitions that fit the same model (BIC-chosen scheme, the successor to PartitionFinder greedy). Pair with `-rcluster 10` (relaxed clustering: only test the top 10% most-similar pairs) for many partitions.
- **Site-heterogeneous mixtures for deep data.** Empirical matrices (LG, WAG) assume one residue-frequency vector for the whole alignment; real proteins do not, and that across-site compositional heterogeneity is the chief driver of deep LBA. The ML answer is the C10..C60 profile-mixture series (`LG+C60+F+G`), made tractable by PMSF (Wang 2018): a guide-tree pass computes one posterior-mean profile per site, and the real search uses those frozen profiles. PMSF is the standard recommendation for deep / LBA-prone protein phylogenomics.

```bash
# PMSF (Wang 2018): a reproducible guide tree, then C-mixture profiles fit on that fixed topology
iqtree3 -s aln.fasta -m LG+F+G --seed 12345 --prefix guide
iqtree3 -s aln.fasta -m LG+C60+F+G -ft guide.treefile -B 1000 -bnni --seed 12345 --prefix pmsf
#  -ft guide.treefile   fixes the topology for the PMSF profile fit; it does not run a second topology search
```
The `.iqtree` report from step 2 shows `Model of substitution: LG+SSF+F+G4` with "site specific frequencies" --
that SSF (site-specific frequency) model, cited to Wang 2018, IS the fitted PMSF profile derived from `LG+C60+F+G`;
IQ-TREE reports the fitted site-specific frequencies rather than echoing the mixture name verbatim. On a
15-taxon/445-site protein alignment C60+UFBoot1000 took ~4 minutes single-threaded (`-T AUTO` for more cores);
drop to `C20`/`C10` if `C60` is too slow for the taxon count and add `-T AUTO`.

## Per-Method Failure Modes

### Long-Branch Attraction
**Trigger:** Two or more independently fast-evolving lineages on long branches separated by a short internode.
**Mechanism:** Convergent/homoplastic substitutions on the long branches look like shared ancestry; a site-homogeneous model cannot separate convergence from homology and groups them -- bias that GROWS with more sites.
**Symptom:** Fast taxa group with the outgroup or each other at 100% bootstrap; the grouping collapses under a better model or when a long-branch taxon is removed.
**Fix:** Site-heterogeneous model (C60/PMSF) first; remove the fastest sites and watch the node; drop the long-branch taxon or use a closer outgroup; cross-check with SR4/Dayhoff recoding. Believe a deep node only when it survives all of these, not when it merely has UFBoot 100 under LG+G.

### Model Underspecification
**Trigger:** A single inadequate model on heterogeneous data; the best site-homogeneous model by BIC still inadequate at depth.
**Mechanism:** Wrong matrix, missing rate heterogeneity, or no partitioning biases the topology while support stays high.
**Symptom:** Biologically implausible nodes with full support that move under a richer model.
**Fix:** `-m MFP` (+MERGE for multi-locus), FreeRate `+R`, and at amino-acid depth a C60/PMSF mixture. Best-by-BIC among site-homogeneous models is not sufficient deep in the tree.

### Over-Partitioning
**Trigger:** Hundreds of hand-defined partitions, especially under `-Q`.
**Mechanism:** Each partition's model is estimated from too little data; parameter and branch-length estimates get noisy without reducing bias.
**Symptom:** Slow runs, noisy estimates, degraded support.
**Fix:** `-m MFP+MERGE -rcluster 10` to the coarsest BIC-justified scheme; use `-p`, not `-Q`.

### The Support-Accuracy Gap
**Trigger:** UFBoot/bootstrap ~100 everywhere, including implausible or conflicting nodes.
**Mechanism:** Support measures repeatability of a possibly-biased estimate; concatenation pools conflicting gene signals so the likelihood is certain while the loci disagree.
**Symptom:** Full support but gCF ~33 / sCF near its ~33% floor on contested nodes.
**Fix:** Compute gCF/sCFL; treat UFBoot 100 + gCF ~33 as UNRESOLVED; require SH-aLRT >=80 AND UFBoot >=95; use `-bnni`. If most genes reject the ML resolution, route to a coalescent species tree -> species-trees.

## Quantitative Thresholds

| Quantity | Threshold | Source |
|----------|-----------|--------|
| UFBoot strong support | >=95 (NOT the bootstrap 70) | Hoang 2018 |
| SH-aLRT strong support | >=80 | Guindon 2010 |
| Joint rule | SH-aLRT >=80 AND UFBoot >=95 | community standard |
| Standard bootstrap "good" | >=70 (different metric, folklore) | Felsenstein 1985 / Hillis 1993 |
| aBayes strong | >=0.95 (anti-conservative; never alone) | Anisimova 2011 |
| Bootstrap replicates | UFBoot `-B` >=1000; SH-aLRT `-alrt` >=1000; AU `-zb` 10000 | IQ-TREE docs |
| gCF reading | >~75 agree; ~50 conflicted; ~33 effective polytomy; <33 with higher alternative = possible wrong resolution | Minh 2020 |
| sCF floor | ~33% (three quartet resolutions); ~33 = no site signal | Minh 2020 |
| AU test | p-AU < 0.05 => tree REJECTED | Shimodaira 2002 |
| Model selection | rank by BIC (`-m MFP` default), k*ln(n) penalty | Kalyaanamoorthy 2017 |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Reading UFBoot 80 as "supported" | applied the bootstrap-70 rule to a different scale | use UFBoot >=95 AND SH-aLRT >=80 |
| Fully-supported deep node distrusted by reviewer | no concordance factors reported | compute gCF/sCFL; treat high-support/low-CF as unresolved |
| `--scfl` unrecognized | IQ-TREE older than 2.2.2 | upgrade, or fall back to parsimony `--scf` |
| "Do not specify --scf or --gcf with --scfl" | gCF and likelihood sCF in one call | run `-t ... --gcf` and `-te ... -s ... --scfl` as two calls |
| `iqtree2: command not found` after `conda install -c bioconda iqtree` | bioconda now ships IQ-TREE 3 | call `iqtree3` (or `iqtree`) |
| Concatenated tree confidently wrong on a rapid radiation | ILS; concatenation is inconsistent in the anomaly zone | infer per-locus gene trees and a coalescent species tree -> species-trees |
| AU test "fails to reject" the alternative | weak data, or candidate set padded | do not pad the set; failure to reject is not acceptance |
