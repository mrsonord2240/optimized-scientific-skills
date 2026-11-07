# Lineage methods and limitations

Source API and numerical results were not independently verified. Source examples use illustrative assumptions; compare solvers and do not infer diagnosis from ancestry.

## Assay Decision Table

Choose the recording technology by the question, not by availability.

| Assay | Mark / model | Use when | Fails when |
|-------|--------------|----------|------------|
| CRISPR scar array (GESTALT/scGESTALT/ScarTrace/LINNAEUS) | cumulative irreversible Cas9 indels = phylogenetic characters | deep tree topology in an engineered organism; scar + scRNA-seq in the same cell | saturation caps depth; homoplasy fuses lineages; heritable dropout deletes clades; not usable in native human tissue |
| Static expressed barcode (LARRY, Weinreb 2020) | one unique inherited lentiviral barcode per founder = a flat clone | clean state->fate maps via split-and-profile; proving state underdetermines fate | gives clonal membership, NOT division-order topology; library << founders causes barcode collisions |
| Combinatorial/sequential barcode (CellTag, Biddy 2018) | combination of expressed tags + nested timepoints | extra clonal resolution and coarse multi-level nesting | shallow trees; collisions; not a resolved phylogeny |
| Somatic mtDNA (Ludwig 2019; mtscATAC-seq; MAESTER) | drifting heteroplasmy of somatic mtDNA variants | retrospective tracing in primary human tissue with no engineering | low mutation rate -> few informative variants; hotspot homoplasy; coverage/dropout; heteroplasmy drift + selection; gives clonal grouping not deep trees |

## Reconstruction Solver Decision Table

For scar data, the solver is a separate choice from the assay (Cassiopeia, Jones 2020; Startle, Sashittal 2023).

| Solver | Model | Use when | Fails when |
|--------|-------|----------|------------|
| VanillaGreedySolver | top-down parsimony, split on most-frequent mutation | fast first pass; 10^4-10^5 cells | greedy split errors propagate; sensitive to homoplasy |
| ILPSolver | integer-LP Steiner tree (Gurobi); near-optimal | small clades needing accuracy | expensive; does not scale to large trees; needs a separately licensed Gurobi install |
| HybridSolver | greedy top + ILP on small subclades | the practical default for large data | inherits greedy errors at the top split |
| NeighborJoiningSolver | distance-based (weighted Hamming) | quick comparison baseline; non-character distances | less accurate than parsimony on scar characters |
| Startle (Startle-ILP / Startle-NNI) | star-homoplasy: a character mutates at most once per root-to-leaf path | severe homoplasy/dropout breaks parsimony | ILP cost; NNI is heuristic at scale |

Run a panel of solvers, not one: it is rare for a single solver to be optimal over all parts of a tree, and agreement across solvers is the practical certainty signal. Always weight indels by their formation probability (down-weight frequent low-information scars) and report robustness to homoplasy and dropout. Methodology evolves; verify the current solver API and recommended defaults against the installed Cassiopeia docs.

### Compare Solvers and Score Tree Robustness

**Goal:** Quantify how much the topology depends on solver choice and on homoplasy/dropout.
**Approach:** Solve with several solvers, then compare the resulting trees with Robinson-Foulds and the depth-stratified triplets-correct metric.

```python
# bottom_solver=ILPSolver() needs a licensed Gurobi install (see Installation); swap in
# VanillaGreedySolver() as the bottom_solver for a license-free greedy/greedy fallback
hybrid = cas.solver.HybridSolver(top_solver=cas.solver.VanillaGreedySolver(), bottom_solver=cas.solver.ILPSolver(), cell_cutoff=200)
# add_root=True is required by the installed API, not optional -- omitting it raises DistanceSolverError
nj = cas.solver.NeighborJoiningSolver(dissimilarity_function=cas.solver.dissimilarity_functions.weighted_hamming_distance, add_root=True)
for s in (hybrid, nj):
    s.solve(tree)                                       # solve independent copies in practice
rf, rf_max = cas.critique.robinson_foulds(tree_a, tree_b)
triplet_acc = cas.critique.triplets_correct(tree_a, tree_b)
```

Triplets-correct is depth-stratified, so it exposes the field's hard truth: deep (near-root) splits are the least certain and the most consequential, while well-supported leaf structure is often the least interesting biologically.

### Build a Character Matrix From Raw Reads

**Goal:** Go from aligned barcode reads to an allele table and character matrix.
**Approach:** Resolve UMIs, align to the reference, call alleles, group cells into clonal populations, then convert the allele table.

```python
umi_table = cas.pp.resolve_umi_sequence(molecule_table, output_directory='.', min_umi_per_cell=10)
aligned = cas.pp.align_sequences(umi_table, ref_filepath='barcode_reference.fa')
alleles = cas.pp.call_alleles(aligned, ref_filepath='barcode_reference.fa')
alleles = cas.pp.call_lineage_groups(alleles, output_directory='.')
char_matrix, priors, state_map = cas.pp.convert_alleletable_to_character_matrix(alleles)
```

`convert_alleletable_to_character_matrix` returns indel priors alongside the matrix; pass those priors to the solver so frequent low-information scars are down-weighted against homoplasy.

### Group Clones from mtDNA Heteroplasmy

**Goal:** Turn a cells x variants heteroplasmy matrix into a clonal grouping — retrospective tracing recovers clonal blobs, never a deep ordered tree (see Assay Decision Table).
**Approach:** Blacklist hotspot/NUMT-like variants that are common regardless of clone, then cluster the remaining heteroplasmy profile and pick the clone count by silhouette score.

```python
import pandas as pd
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import pdist
from sklearn.metrics import silhouette_score

het = pd.read_csv('mtdna_heteroplasmy.csv', index_col=0)   # cells x variants, values in [0, 1]

# A variant present (heteroplasmy > 0.1) in most cells cannot discriminate clones --
# treat it as a hotspot/NUMT signature and drop it rather than trusting it as a mark.
prevalence = (het > 0.1).mean(axis=0)
het_f = het.loc[:, prevalence <= 0.8]

dist = pdist(het_f.values, metric='correlation')
Z = linkage(dist, method='average')

# Choose the clone count by silhouette score over a small candidate range
best_k, best_score, best_labels = None, -1, None
for k in range(2, 8):
    labels = fcluster(Z, t=k, criterion='maxclust')
    if len(set(labels)) < 2:
        continue
    score = silhouette_score(het_f.values, labels, metric='correlation')
    if score > best_score:
        best_k, best_score, best_labels = k, score, labels

clones = pd.Series(best_labels, index=het.index, name='mtdna_clone')
```

Verified on a synthetic 40-cell/10-variant heteroplasmy matrix with 4 planted clones plus
one recurrent hotspot variant: the blacklist step drops exactly the hotspot, and
clustering recovers the 4 planted clones exactly (adjusted Rand index 1.0). Report this
as a clonal grouping only — never branch order or divergence times. mtDNA heteroplasmy from primary human tissue is human genomic data: handle it under the same consent and privacy requirements as any other.

### Refine a Tree Under Severe Homoplasy with Startle

**Goal:** Re-optimize a seed tree's topology under the star-homoplasy model (each character mutates at most once per root-to-leaf path), for data too homoplastic for greedy parsimony to trust.
**Approach:** Build a fast seed tree (`VanillaGreedySolver`), export it and the character matrix, and let Startle's NNI hill-climber (CLI, not a Python import) search for lower-weighted-parsimony topologies.

```python
solver = cas.solver.VanillaGreedySolver()
solver.solve(tree, collapse_mutationless_edges=True)
with open('seed_tree.newick', 'w') as f:
    f.write(tree.get_newick())

# Startle wants per-character, per-state mutation priors as a long CSV. From raw reads,
# convert_alleletable_to_character_matrix already returns these (see above); from a
# plain character matrix, estimate them empirically:
rows = []
for c in char_matrix.columns:
    col = char_matrix[c]
    edited = col[(col != 0) & (col != -1)]
    for state, p in edited.value_counts(normalize=True).items():
        rows.append({'character': c, 'state': int(state), 'probability': float(p)})
pd.DataFrame(rows).to_csv('startle_priors.csv', index=False)
char_matrix.to_csv('startle_char_matrix.csv')   # row index must match the Newick leaf names
```

```bash
startle large startle_char_matrix.csv startle_priors.csv seed_tree.newick \
    --output refined --iterations 100
```

`refined_tree.newick` is the re-optimized topology; `refined_info.json` reports the
weighted parsimony score per iteration, so you can confirm the search actually improved
on the seed tree instead of trusting the exit code. Verified on this Skill's own
60-cell/16-site synthetic scar matrix: a VanillaGreedySolver seed tree at weighted
parsimony ~7.5 was refined to 6.59 over 20 NNI iterations.

## Threshold and Parameter Rationale

| Parameter | Typical value | Rationale |
|-----------|---------------|-----------|
| min_umi_per_cell | ~10 | below this, allele calls are dominated by sequencing noise |
| missing fraction per cell | drop > ~0.5 | cells missing most characters carry little phylogenetic signal and inflate ambiguity |
| informative character | states in > 1 cell | a scar seen in one cell cannot group lineages; uninformative for topology |
| indel prior weighting | from empirical indel frequencies | frequent microhomology-driven indels are high-homoplasy, low-information; down-weight them |
| barcode library complexity | >> number of founders | small libraries cause collisions (two founders share a barcode -> phantom merged clone) |
| HybridSolver cell_cutoff | ~200 | subclades below the cutoff are solved exactly by ILP; above it, greedily |

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| Distinct lineages collapse into one clade | missing data coded as the unedited 0 state | keep -1 missing distinct from 0; model dropout, never treat it as unedited |
| Parsimony fuses unrelated cells | homoplasy: independent cells share a frequent indel | weight indels by formation probability; use Startle's star-homoplasy model under heavy convergence |
| Late divisions are unresolved near the leaves | editable array saturated; recording stopped early | use inducible/paced recorders; report the per-site edit fraction distribution |
| Two founders appear as one giant clone | barcode library too small relative to founders -> collision | use library complexity >> cell number; estimate collisions empirically |
| Impossible chimeric clones / character vectors | doublets carry two barcode/scar sets | run doublet detection and barcode-consistency filtering before reconstruction |
| Deep splits flip between solvers | early splits rest on the fewest, most-overwritten characters | report branch support; trust leaf structure more than the root; run a solver panel |
| mtDNA "tree" is actually clonal blobs | low somatic mutation rate; hotspot homoplasy; heteroplasmy drift and selection | claim clonal grouping not deep ordered trees; blacklist NUMTs/RNA-edit/hotspot sites |
| State-based branch call confidently wrong | state underdetermines fate (Weinreb 2020); map is one-to-many | frame fate as a prediction; validate with prospective lineage data, integrate with CoSpar |
| ILPSolver/HybridSolver fails with a cryptic, unrelated-looking error | `gurobipy` not installed or licensed | use `HybridSolver(bottom_solver=VanillaGreedySolver())`, or license Gurobi |
| Tree/clone output used to claim malignancy or treatment need | a clonal/lineage tree shows ancestry and clonal relationships only, not histopathology | decline the diagnostic or treatment claim; recommend pathologist/oncologist review |

## References

Weinreb C, Rodriguez-Fraticelli A, Camargo FD, Klein AM (2020). Lineage tracing on transcriptional landscapes links state to fate during differentiation (LARRY). Science 367(6479):eaaw3381.
McKenna A, Findlay GM, Gagnon JA, Horwitz MS, Schier AF, Shendure J (2016). Whole-organism lineage tracing by combinatorial and cumulative genome editing (GESTALT). Science 353(6298):aaf7907.
Raj B, Wagner DE, McKenna A, et al. (2018). Simultaneous single-cell profiling of lineages and cell types in the vertebrate brain (scGESTALT). Nat Biotechnol 36(5):442-450.
Alemany A, Florescu M, Baron CS, Peterson-Maduro J, van Oudenaarden A (2018). Whole-organism clone tracing using single-cell sequencing (ScarTrace). Nature 556(7699):108-112.
Spanjaard B, Hu B, Mitic N, et al. (2018). Simultaneous lineage tracing and cell-type identification using CRISPR-Cas9-induced genetic scars (LINNAEUS). Nat Biotechnol 36:469-473.
Biddy BA, Kong W, Kamimoto K, et al. (2018). Single-cell mapping of lineage and identity in direct reprogramming (CellTag). Nature 564:219-224.
Wang SW, Herriges MJ, Hurley K, Kotton DN, Klein AM (2022). CoSpar identifies early cell fate biases from single-cell transcriptomic and lineage information. Nat Biotechnol 40:1066-1074.
Ludwig LS, Lareau CA, Ulirsch JC, et al. (2019). Lineage tracing in humans enabled by mitochondrial mutations and single-cell genomics. Cell 176(6):1325-1339.
Lareau CA, Ludwig LS, Muus C, et al. (2021). Massively parallel single-cell mitochondrial DNA genotyping and chromatin profiling (mtscATAC-seq). Nat Biotechnol 39:451-461.
Miller TE, Lareau CA, Verga JA, et al. (2022). Mitochondrial variant enrichment from high-throughput single-cell RNA sequencing resolves clonal populations (MAESTER). Nat Biotechnol 40:1030-1034.
Jones MG, Khodaverdian A, Quinn JJ, et al. (2020). Inference of single-cell phylogenies from lineage tracing data using Cassiopeia. Genome Biology 21:92.
Sashittal P, Schmidt H, Chan M, Raphael BJ (2023). Startle: a star homoplasy approach for CRISPR-Cas9 lineage tracing. Cell Systems 14(12):1113-1121.


### Integrate Clones With State Using CoSpar

**Goal:** Recover early fate bias from sparse clonal barcodes rather than assuming the manifold encodes fate.
**Approach:** Fit a transition map jointly from clonal observations and transcriptomic similarity, then read fate bias and fate maps.

```python
import cospar as cs, scanpy as sc
adata = cs.hf.read('lineage_traced.h5ad')   # X = raw counts; obsm['X_clone'] = cells x clones; obs['time_info']
# CoSpar needs a precomputed PCA and 2-D embedding. initialize_adata_object only WARNS when they are
# missing, and the crash comes later as KeyError: 'X_emb'. Compute them on a log-normalized copy so
# adata.X stays raw counts (CoSpar warns on log-transformed X).
emb = adata.copy()
sc.pp.normalize_total(emb); sc.pp.log1p(emb); sc.pp.pca(emb); sc.pp.neighbors(emb); sc.tl.umap(emb)
# data_des must be unique per dataset: CoSpar caches similarity matrices on disk under it
adata = cs.pp.initialize_adata_object(adata, X_clone=adata.obsm['X_clone'], time_info=adata.obs['time_info'],
                                      X_pca=emb.obsm['X_pca'], X_emb=emb.obsm['X_umap'], data_des='my_dataset')
adata = cs.tmap.infer_Tmap_from_multitime_clones(adata, smooth_array=[15, 10, 5], sparsity_threshold=0.1)
cs.tl.fate_bias(adata, selected_fates=['Monocyte', 'Neutrophil'])
cs.pl.fate_bias(adata, selected_fates=['Monocyte', 'Neutrophil'])
```

`infer_Tmap_from_multitime_clones` was run here with 2 and with 3 time points in `time_info`. If a second dataset raises `ValueError: The pre-computed similarity matrix does not have the right dimension`, it reused the cached matrices of an earlier dataset with the same `data_des` (the default is `'cospar'`). Give each dataset its own `data_des`, or pass `compute_new=True` to the inference call.

CoSpar operationalizes Weinreb 2020: it propagates fate probabilities onto cells lacking clonal labels and is robust to severe downsampling of lineage data, but it needs paired clone + state and does NOT build a phylogenetic tree (clones are flat). CoSpar needs MULTIPLE independent clones to be lineage-informed; with effectively one clone the constraint is vacuous and the transition map degenerates to transcriptomic similarity, the state-only answer CoSpar exists to correct. For tree topology from scars, use Cassiopeia or Startle.

`infer_Tmap_from_multitime_clones` took ~15 s on a 200-cell toy dataset with the default
`smooth_array=[15, 10, 5]` (2026-09-21); cost grows with cell number, so expect much longer on real
datasets and use a smaller `smooth_array` for a first exploratory pass.
