## Branch Support -- the Heart

Five measures, three different questions; the cardinal sin is cross-comparing their cutoffs.

| Metric | What it perturbs / measures | Strong cutoff | Tool / flag | Failure mode |
|--------|------------------------------|---------------|-------------|--------------|
| Standard bootstrap (FBP) | resample sites; clade frequency (binary) | >=70 (folklore) | RAxML-NG `--bs-metric fbp`; IQ-TREE `-b` | slow; crushed by rogue taxa in big trees |
| Ultrafast bootstrap 2 (UFBoot) | RELL-resampled log-Ls; ~unbiased clade prob | **>=95 (NOT 70)** | IQ-TREE2 `-B 1000` (+`-bnni`) | inflates under model violation -> use `-bnni` |
| SH-aLRT | local NNI likelihood ratio (no resampling) | >=80 | IQ-TREE2 `-alrt 1000` | conservative on very short branches |
| aBayes | posterior from 3 NNI Ls, flat prior | >=0.95 | IQ-TREE2 `-abayes` | anti-conservative; never the sole criterion |
| Transfer bootstrap (TBE) | gradual transfer distance under resampling | no fixed cutoff; > FBP | RAxML-NG `--bs-metric tbe` | permissive; "fuzzy" branch identity |

UFBoot2 (Hoang 2018) uses the RELL trick (resample site log-likelihoods, reuse a candidate tree set) to run hundreds of times faster than the Felsenstein bootstrap (1985), and its values are CLOSER to unbiased clade probabilities -- which is exactly why the strong-support cutoff is 95, not the conservative-bootstrap 70. Treating UFBoot 70 as "good" is a category error. `-bnni` re-optimizes each replicate tree by NNI to rein in the inflation that model violation causes; use `-B 1000 -bnni` routinely. UFBoot is on a DIFFERENT scale from the standard bootstrap -- never read it with the BP-70 rule.

SH-aLRT (Guindon 2010) does not resample data; for each branch it tests whether the ML likelihood beats its two best NNI rearrangements. UFBoot (data perturbation) and SH-aLRT (tree perturbation) have different failure modes, so the community-standard joint criterion requires both:

> A branch is strongly supported iff SH-aLRT >= 80% AND UFBoot >= 95%.

For large rogue-taxon-prone trees, the binary Felsenstein bootstrap lets a single wandering tip crush an otherwise-recovered deep branch; transfer bootstrap (TBE, Lemoine 2018, RAxML-NG `--bs-metric tbe`) replaces the in/out indicator with a gradual transfer distance and rescues those branches, at the cost of being more permissive.

## Pre-flight: Small or Low-Signal Alignments

**Goal:** Decide whether an alignment carries enough signal to trust the tree, and whether apparent structure is really a duplicate sequence, before interpreting support. This is the opposite failure mode from the rest of this Skill -- few taxa, short or conserved alignments, not phylogenomic scale.

**Approach:** IQ-TREE2 already reports every diagnostic below in its own `.iqtree`/`.log` output -- do not compute them separately.

- **Parsimony-informative sites.** The `.iqtree` report's SEQUENCE ALIGNMENT block prints `Number of parsimony informative sites: N`. Few informative sites (short or conserved data, few taxa) means little of the alignment actually drove the topology; check this count before trusting a deep or well-supported-looking node. Example, 6 short primate isolates / 250 sites: `Number of constant sites: 187 (= 74.8% of all sites)` / `Number of parsimony informative sites: 33`.
- **Identical sequences are KEPT, not dropped.** IQ-TREE2 does not collapse or remove identical sequences by default; it flags the pair in the log and keeps both in the tree: `NOTE: <seq> is identical to <seq> but kept for subsequent analysis`. The duplicate then sits on its own near-zero branch next to its twin (see next point) -- do not mistake the pair for two independently resolved tips.
- **Collapse or flag branches <=1e-6 before reading support.** IQ-TREE2's branch-length floor is 1e-6 (printed as e.g. `0.0000010000`); a branch at or near that floor means the two sides are effectively identical (a duplicate or near-duplicate), not a resolved split. Support on the node immediately adjacent to a near-zero branch is answering "is this near-duplicate reliably near-duplicate," not "is this clade real" -- collapse such branches first (e.g. `ape::di2multi(tree, tol=1e-6)` in R) before reading support on the adjacent node.

Verified on a 6-taxon alignment with one exact duplicate (`Homo_sapiens` / `Homo_sapiens_isolate2`, `iqtree2 -s isolates6.fa -m MFP -B 1000 -bnni --alrt 1000 -T 1`): `.iqtree` reports `Number of parsimony informative sites: 33` (of 250); `.log` prints `NOTE: Homo_sapiens_isolate2 is identical to Homo_sapiens but kept for subsequent analysis`; the treefile resolves both duplicate tips at `0.0000010000` (IQ-TREE's 1e-6 floor), and the branch grouping that pair with the rest of the tree carries an 86.7/100 SH-aLRT/UFBoot split driven entirely by the duplicate, not by independent phylogenetic signal.

## Concordance Factors

Bootstrap quantifies statistical confidence given the concatenated data; concordance factors (Minh 2020) quantify how much of the actual data carries a branch.

- **gCF (gene concordance factor):** the percentage of decisive single-locus gene trees that contain the exact branch. Needs per-locus gene trees.
- **sCF (site concordance factor):** the percentage of decisive sites supporting the branch, from sampled quartets; works on a single concatenated alignment with no gene trees.
- **sCFL (likelihood sCF, Mo 2023):** uses ancestral-state likelihoods rather than parsimony quartet counting, substantially reducing (not abolishing) homoplasy and taxon-sampling bias. Prefer `--scfl` over the old `--scf` on IQ-TREE 2.2.2+.

```bash
# one gene tree per locus from a directory of locus alignments (-S = separate, no concatenation)
iqtree3 -S loci_dir -m MFP -B 1000 -T AUTO --prefix loci   # -B 1000 = UFBoot per gene tree, needed to contract weak branches before ASTRAL

# gene + likelihood site concordance against the fixed concatenated tree: TWO calls,
# IQ-TREE rejects --gcf together with --scfl ("Do not specify --scf or --gcf with --scfl")
iqtree3 -t concat.treefile --gcf loci.treefile --prefix concord_g                 # per-locus gene trees for gCF
iqtree3 -te concat.treefile -s concat.fasta --scfl 100 -T 4 --prefix concord_s    # 100 quartets per branch for likelihood sCF
```

Outputs `concord_g.cf.tree` / `concord_g.cf.stat` (gCF labels; per-branch gCF, gDF1, gDF2, gDFP) and `concord_s.cf.tree` / `concord_s.cf.stat` (likelihood sCF). A node with UFBoot 100 but gCF ~35 (gDF1 ~33, gDF2 ~30) is genes split three ways: the concatenated point estimate barely edges the alternatives and the node is biologically unresolved -- the signature of ILS or introgression, not a clade. Report CFs alongside support on every phylogenomic tree.

## Topology Tests

For testing an a-priori hypothesis ("can I reject that X and Y are monophyletic?") against the ML tree, by comparing a set of fixed trees. Build the constrained tree with `-g constraint.tree`, then evaluate both trees:

```bash
# trees.nex holds the unconstrained ML tree + the constrained/alternative trees
iqtree3 -s aln.fasta -m <model> -z trees.nex -n 0 -zb 10000 -au --prefix autest
#  -z trees.nex   trees to compare        -n 0   no fresh search, just evaluate
#  -zb 10000      RELL replicates (>=1000) -au    add the AU test (must accompany -zb)
```

The AU test (Shimodaira 2002) uses multiscale bootstrap resampling to correct both the selection bias of the KH test (invalid on the data-selected ML tree) and the over-conservatism of the SH test (which rejects less as the candidate set is padded). Use AU by default. **p-AU < 0.05 means that tree is REJECTED**; p-AU >= 0.05 means it is in the 95% confidence set (failure to reject is not acceptance -- weak data fails to reject many trees). Report SH/KH only for completeness.

## Partitioning

When splitting an alignment into partitions, the branch-length linkage choice is what people get wrong:

| Mode | Flag | Branch lengths | Use when |
|------|------|----------------|----------|
| Edge-equal | `-q part.nex` | identical across partitions | partitions share rate (rare, restrictive) |
| Edge-linked proportional | `-p part.nex` | shared topology, per-partition rate multiplier | DEFAULT -- genes evolve at different speeds, share history |
| Edge-unlinked | `-Q part.nex` | fully independent per partition | genuine heterotachy; parameter-hungry, overfits |

`-p` (edge-linked proportional) is the standard: one rate multiplier per partition over a shared topology. Over-partitioning spends degrees of freedom without bias reduction and inflates variance; the antidote is to start fine (gene x codon position) and let `-m MFP+MERGE -rcluster 10` find the coarsest BIC-justified scheme. Prefer a merged scheme over a hand-picked maximal one; prefer `-p` over `-Q`.
