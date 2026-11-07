## Per-Method Failure Modes

### Distant Outgroup Roots Inside the Ingroup
**Trigger:** A single distant outgroup, or many outgroups all far from the ingroup, used to root.
**Mechanism:** The long outgroup branch attracts the fastest ingroup taxon (LBA), pulling the root into the ingroup; far outgroups attach at essentially random positions (DeSalle et al. 2023).
**Symptom:** Ingroup not recovered as monophyletic; a fast taxon appears "basal"; the root jumps as the outgroup set changes.
**Fix:** Use multiple CLOSELY-related, monophyletic, roughly-equidistant outgroups; verify ingroup monophyly; treat a lonely distant outgroup as a red flag.

### Single Distant Outgroup -- Long-Branch Misrooting
**Trigger:** Rooting on one distant outgroup taxon.
**Mechanism:** One unbroken long branch maximizes LBA exposure, with no way to subdivide it or check ingroup monophyly, so the root is drawn toward other long branches.
**Symptom:** The root lands inside the ingroup or on a spurious deep branch; deep nodes near the root are unstable across analyses.
**Fix:** Add multiple closer, monophyletic outgroups to subdivide the long branch and enable a monophyly check, or use an outgroup-free method (MAD/MinVar).

### Midpoint Misroots Under Rate Variation
**Trigger:** Midpoint rooting a deep tree with heterotachy or any long branch.
**Mechanism:** Midpoint assumes a clock; a long branch hijacks the longest path and slides the root onto the fast lineage.
**Symptom:** A rate-elevated taxon appears earliest-diverging; the root sits on a suspiciously long branch.
**Fix:** Use MAD or MinVar (clock-relaxed) or a close outgroup; report root uncertainty and hedge "basal".

### Soft vs Hard Polytomy Conflation
**Trigger:** Collapsing branches below a support threshold, then describing the multifurcation as a radiation.
**Mechanism:** Threshold-collapse encodes UNCERTAINTY (soft polytomy = "we cannot resolve the order"); a hard polytomy is a biological claim of simultaneous divergence -- different meanings.
**Symptom:** A paper claims simultaneous divergence from what is just unresolved data.
**Fix:** Label threshold-collapsed nodes as unresolved/soft; never read them as a biological radiation.

### Pruning Leaves a Spurious Node or Drops Distances
**Trigger:** A tool that leaves the degree-2 node in place, or removes it without summing branch lengths (notably ete3 `prune` without `preserve_branch_length=True`).
**Mechanism:** The suppressed knee's branch length is not added to its child, so every path through that lineage shortens.
**Symptom:** Patristic distances shrink; subsequent midpoint/MAD/MinVar rooting on the pruned tree is now wrong.
**Fix:** Use Bio.Phylo `prune` / ape `drop.tip` (correct by default) or pass ete3 `preserve_branch_length=True`; verify a known pairwise distance is unchanged.

### Randomly Resolving a Polytomy Biases Downstream
**Trigger:** `multi2di` / `resolve_polytomy` to satisfy a binary-tree requirement, then analyzing the single tree.
**Mechanism:** An arbitrary order with zero-length branches is invented; the downstream result depends on a topology the data never supported.
**Symptom:** Results that change under a different random resolution.
**Fix:** Integrate over many resolutions and summarize; treat zero-length branches as "no info", not instantaneous divergence.

### Ladderizing Manufactures an Apparent Trend
**Trigger:** Ladderizing a figure so one lineage sits visually at the top/bottom.
**Mechanism:** Rotation about internal nodes changes no topology, branch lengths, or bipartitions, but readers unconsciously read the staircase as an early-to-late sequence.
**Symptom:** Reviewers infer a basal-to-derived trend that does not exist.
**Fix:** Ladderize only for legibility; state that rotation changes no biology. -> tree-visualization.

## Quantitative Thresholds

| Quantity | Value | Source / rationale |
|----------|-------|--------------------|
| Midpoint root recovery, single-outgroup source data | ~54% (barely a coin flip) | Hess and Russo 2007 |
| Midpoint root recovery, multiple-outgroup source data | ~82% (inconsistent) to ~94% (consistent) | Hess and Russo 2007 |
| MAD accuracy on benchmarks | >~70%, beating midpoint | Tria et al. 2017 |
| MinVar vs midpoint | matched or beat midpoint in all simulated conditions | Mai et al. 2017 |
| Bootstrap collapse cutoff | <50% (near-uninformative floor) or <70% (reliability boundary); state which | common practice |
| UFBoot2 collapse cutoff | <95% (a different scale from bootstrap -- not 70) | Hoang 2018 |
| Bayesian posterior collapse cutoff | <0.95 | common practice |
| Zero-length collapse tolerance | match the program that wrote the tree: IQ-TREE writes zero as 1e-6 (use `<= 1e-6`); ape `di2multi` default 1e-8 | a tolerance below the writer's floor catches no zero-length branch |
| Outgroup distance | accuracy degrades monotonically with distance | prefer the closest credible outgroup |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| `ape::midpoint` function-not-found | midpoint is not in base ape | use `phangorn::midpoint` or `phytools::midpoint.root` |
| Pruned tree has wrong distances | degree-2 node not suppressed / not summed | Bio.Phylo `prune` / ape `drop.tip` default, or ete3 `preserve_branch_length=True` |
| "Extract the clade of X,Y,Z" returns extra taxa or None | X,Y,Z are not monophyletic | prune to the taxon set (induced subtree) instead of extracting a clade |
| `di2multi(tol)` did not drop low-support nodes | `tol` filters branch LENGTH, not support | zero out low-support branch lengths first, then `di2multi` |
| Bootstrap re-read as support for a pruned subtree | support was computed on the full taxon set | re-run inference on the subset for valid support |
| MAD and MinVar disagree on the root | the root is genuinely poorly determined | treat as uncertain; seek a close outgroup or a RootDigger (`rd`) confidence run; hedge all ancestral claims |
| `root_with_midpoint` AttributeError | no such method | the methods are `root_at_midpoint()` and `root_with_outgroup()` |
| Valid outgroup reported "not monophyletic" on an IQ-TREE tree | `is_monophyletic` tests rooted clades; the file is rooted on the outgroup's own trifurcation | root on an ingroup tip first, then test and root with `outgroup_branch_length` |
| `collapse_all` on support collapses nothing, no error | IQ-TREE `SH-aLRT/UFBoot` label kept in `clade.name`, `confidence` None | parse `clade.name.split('/')[-1]`; assert support was read before collapsing |
