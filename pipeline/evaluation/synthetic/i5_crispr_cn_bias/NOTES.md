# i5_crispr_cn_bias (synthetic; not agent-visible)

## Pitfall and source
Cas9 cuts at every copy of an amplified locus; the resulting DNA-damage toxicity depletes guides against
non-essential genes in high-copy amplicons, so they look like dependencies unless copy number is accounted for
(Aguirre et al. 2016 Cancer Discov; Meyers et al. 2017 Nat Genet, which also proposed CN-correction methods).

## Construction (`generate.py`, seed 20260919, deterministic, ~1 s)
- 600 invented genes (SYN####, names shuffled w.r.t. position and role) x 4 guides + 100 NTC guides; 12 synthetic chromosomes.
- Per-guide true log2FC `L = e * (E_gene + C_cn)`, `e ~ Beta(4,1.5)` (efficiency), plasmid lognormal skew (sigma 0.5),
  NB counts (plasmid ~600x, late replicates ~500x / 450x, extra per-replicate per-guide noise sd 0.2 log2).
- Truth = 40 strong essentials (`E` ~ -U(2.4, 5.0)); 2 of them sit inside the amplicon (CN 10.4), so a naive
  exclude-amplicon approach loses 2 (plus any in other gains; 6 essentials have CN > 2.5) but a CN-model keeps them.
- 20 weak genes (`E` ~ -U(0.3, 0.6), gene-level ~ -0.3 to -0.4): below the -1 floor stated in the task, so not truth.
- Amplicon: 32 contiguous chr7 genes in 10 Mb (30 non-essential + 2 essential); segment CN 5.5 / 6.9 / 10.4 / 8.6 / 5.7
  (mean ~7.4). `C = -1.2 * log2(CN/2)` -> gene-level raw LFC about -1.1 (CN 5.5) to -2.2 (CN 10.4), i.e. comparable to real essentials
  (essentials ~ -1.5 to -3.7). Other gains: chr3 arm CN 3 (~-0.4 raw, below the floor), chr10 focal CN 4 (8 genes, ~-0.9 raw,
  borderline). Losses: chr5 arm CN 1 and chr12 focal CN 1, `C = +0.3*log2(2/CN)`.
- Data has no essentiality info; NTC labelled `NTC` in library.tsv. gene_coordinates has no NTC row (must be handled).

## Answer convention (in task text)
Gene-level log2FC (late vs plasmid) about -1 or stronger, 5% FDR, "for reasons intrinsic to loss of that gene rather than
an artefact of the experiment". The -1 floor is what makes weak genes non-truth; it is stated so the truth is not ambiguous.

## Grader
`marker_gene_precision_recall`, flat list, k = len(submitted list) (no `expected_count`), precision >= 0.75, recall >= 0.75.
Case-insensitive; duplicates collapse in the TP set but count in k.

## Results (`python solvers.py`)
| solver | n | TP | precision | recall | pass |
|---|---|---|---|---|---|
| correct: robust regression LFC~log2CN (gain+loss terms) | 40 | 40 | 1.00 | 1.00 | yes |
| correct: CRISPRcleanR-like segment-median subtraction | 40 | 40 | 1.00 | 1.00 | yes |
| correct: exclude genes with CN >= 2.5 | 34 | 34 | 1.00 | 0.85 | yes |
| trap: naive gene LFC <= -1 + FDR (28 amplicon FPs) | 72 | 40 | 0.56 | 1.00 | no |
| trap: RRA-like FDR only, no floor | 100 | 40 | 0.40 | 1.00 | no |
| trap: top-N by LFC, N = #FDR hits | 100 | 40 | 0.40 | 1.00 | no |
| info: CN-corrected but ignoring the floor (adds 5 weak genes) | 46 | 40 | 0.87 | 1.00 | yes |

Ad-hoc variants also pass (both 1.00/1.00): single-slope OLS of guide LFC on log2CN; gene-level gains-only fit; mean/z-test caller.
"Exclude CN >= 4 only" gives 0.90/0.95 (passes). Naive mean/z-test caller: 0.57 (fails).
Regression recovers slope -0.88 gain (design: -1.2 x mean efficiency 0.73 = -0.88) and +0.25 loss shift (design 0.3 x 0.73 = 0.22).

## Tolerance rationale
Thresholds 0.75/0.75 sit midway between the worst legitimate variant (recall 0.85 for exclusion; precision 0.87 for no-floor)
and the best trap (precision 0.56 naive-with-floor). Margins: correct >= 0.10 (recall of exclusion) to 0.25; trap 0.19 below.

## Blind spots / ambiguity / gaming
- BLIND SPOT (recorded in solvers.py, expected to pass): naive gene LFC with a much stricter floor of -2 (37 genes, precision 0.86,
  recall 0.80) passes; so does a naive LFC ranking truncated at an arbitrary N of about 35-45 (top-40: 0.80/0.80). Both need an
  effect-size cut the task does not state (it says about -1) or luck on N; at the stated -1 floor naive gives 0.56, top-N-by-FDR gives 0.40.
  Raising precision above 0.86 would fail legitimate CN-corrected variants with a few weak-gene FPs, so thresholds stay at 0.75.
- The task tells the solver the floor is -1; a solver that skips the floor after CN correction still passes (0.87), by design.
- Without the FDR/floor both stated in the task, calls differ; convention is explicit, so ambiguity is low.
- Guide efficiency is not observable; sequences are random and uninformative.
- The two amplicon essentials are found only by CN-aware methods that retain them; exclusion-based ones lose them (recall 0.85, still passes).
- An agent could recognise the amplicon from LFC vs genomic position alone even without the CN file (positional smoothing); that is a
  legitimate CN-free route and would also pass.
- Skills read only after the design was frozen: `relevant_skills` lists the CN-correction, MAGeCK, hit-calling, BAGEL, QC and pipeline skills.
