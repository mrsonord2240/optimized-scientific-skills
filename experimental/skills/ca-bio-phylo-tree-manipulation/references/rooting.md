## Rooting Method Decision

| Method | Core assumption | Needs | Robust to rate variation | Gives root uncertainty | Best when | Avoid when |
|--------|-----------------|-------|--------------------------|------------------------|-----------|------------|
| Outgroup (multiple, close) | outgroup truly outside ingroup, not LBA-misplaced | a priori outgroup taxa + sequences | moderate (depends on outgroup branch) | weakly (via ingroup-monophyly check) | closely related, balanced, monophyletic outgroups exist | only a distant/lonely outgroup available |
| Outgroup (single) | same, but one long branch | one outgroup taxon | poor (max LBA exposure) | no | a close single sister exists, nothing better | only a distant/lonely outgroup (one long branch, max LBA exposure) |
| Midpoint | strict molecular clock | tree with branch lengths | no | no | shallow, clock-like, intraspecific/viral data | deep trees, heterotachy, any long branch |
| MinVar | clock deviations are random/unbiased | tree with branch lengths | better than midpoint | no | outgroup-free, noise looks like unbiased clock scatter | strong lineage-specific rate shifts |
| MAD | minimal relative ancestor deviation | tree with branch lengths | yes (tolerates heterotachy) | yes (root ambiguity index) | outgroup-free deep/prokaryotic trees; default no-outgroup choice | strongly structured rate variation |
| Non-reversible likelihood (RootDigger / IQ-TREE) | non-reversible model carries directional signal | tree + ALIGNMENT; enough data | yes (model-based) | yes (per-branch likelihood confidence) | the alignment is available and a root confidence is needed | little data; weak signal; compute-limited |
| Relaxed clock (BEAST) -> divergence-dating | explicit clock + tree prior | alignment + calibrations/tip dates; MCMC | yes | yes (posterior over roots) | time-scaled / dated / phylodynamic analyses | a quick structural edit, no dating intended |

Rule: never report a root from a single method without a sanity check. If a close, monophyletic outgroup is available, use it. If not, run MAD and MinVar; disagreement, or a MAD ambiguity index near 1, means the root is poorly determined and all "basal" claims must be hedged. Agreement is NOT root confidence: under a strong clade-wide rate shift both can agree on the same wrong root. With an alignment and a need for a confidence value, use RootDigger (binary `rd`; Bettisworth and Stamatakis 2021) or non-reversible IQ-TREE; commands for both are in the midpoint-fallback recipe below. MAD = Tria et al. 2017; MinVar = Mai et al. 2017; outgroup-free CLIs and Newick Utilities = Junier and Zdobnov 2010.

## Root with an Outgroup (Multiple, Monophyletic)

**Goal:** Root the tree using known sister-group taxa, preferring multiple close outgroups and verifying ingroup monophyly first.

**Approach:** Confirm the outgroup taxa form a monophyletic group, root on the branch separating them from the ingroup, then check a SEPARATELY-declared, a-priori ingroup taxon list is recovered as monophyletic -- if not, the rooting is suspect. That ingroup list must be the researcher's actual focal taxa, not "every tip not currently called outgroup": `root_with_outgroup` always splits the tree into exactly two clades at the new root, so a complement-of-outgroup check is trivially monophyletic by construction and catches nothing -- it does not detect a long-branch-attracted taxon (e.g. a fast-evolving ingroup member) wrongly named as part of the outgroup. An inferred tree is arbitrarily rooted (IQ-TREE writes `(OutA,OutB,ingroup...)`), so first root on an ingroup tip, or the monophyly test returns False for a valid outgroup; pass `outgroup_branch_length`, or Biopython makes the outgroup MRCA a trifurcating root.

```python
from Bio import Phylo

tree = Phylo.read('tree.nwk', 'newick')
outgroup = [{'name': 'OutA'}, {'name': 'OutB'}]      # multiple close outgroups beat a single long branch
ingroup = [{'name': 'I1'}, {'name': 'I2'}, {'name': 'I3'}]   # the actual focal taxa, declared a priori -- NOT "every tip not in outgroup"

tree.root_with_outgroup({'name': 'I1'})              # any INGROUP tip: makes the outgroup a clade of the temporary root
if tree.is_monophyletic([tree.find_any(**o) for o in outgroup]):
    stem = tree.common_ancestor(*outgroup).branch_length
    tree.root_with_outgroup(*outgroup, outgroup_branch_length=stem / 2)   # bifurcating root on the outgroup stem
else:
    print('outgroup not monophyletic: root placement is unreliable, re-check taxon choice')

if not tree.is_monophyletic([tree.find_any(**i) for i in ingroup]):
    print('ingroup not monophyletic after rooting: outgroup choice is likely wrong (e.g. long-branch attraction pulled a fast ingroup taxon into the outgroup)')
```

## Root at the Midpoint (Clock-Limited Fallback)

**Goal:** Root an outgroup-free tree where evolution is approximately clock-like (shallow/viral data).

**Approach:** Place the root at the midpoint of the longest tip-to-tip path; trust it only when no single long branch can hijack that path. For deep trees with rate variation prefer MAD/MinVar.

```python
tree = Phylo.read('tree.nwk', 'newick')
tree.root_at_midpoint()                              # assumes a clock; a long branch slides the root onto the fast lineage
# Outgroup-free and clock-relaxed (deep trees): standalone CLIs run on the Newick file
# MAD:    mad tree.nwk           -> tree.nwk.rooted   (report the ambiguity index AI; AI near 1 = ambiguous root)
#         a Newick written by Bio.Phylo ends in '):0;', which MAD 2.2 rejects as corrupt -> strip it first: re.sub(r'\):[0-9.eE+-]+;$', ');', nwk)
# MinVar: FastRoot.py -i tree.nwk -m MV -o rooted.nwk
# Non-reversible likelihood root with confidence (needs the alignment; iqtree2 on IQ-TREE 2.x):
#   iqtree3 -s aln.fa --model-joint 12.12 -B 1000 --root-test -zb 1000 -au
#   -> report the rootstrap on the root branch (.rootstrap.nex) and the AU confidence set of root branches (.roottest.csv)
# RootDigger (needs the alignment too; binary is `rd`, not `rootdigger` -- not installed on this machine,
# command checked against the computations/root_digger README and Bettisworth & Stamatakis 2021, not run here):
#   rd --msa aln.fa --tree tree.nwk --exhaustive
#   -> reports the likelihood weight ratio (LWR) per branch; add --early-stop to shorten the exhaustive search
```
