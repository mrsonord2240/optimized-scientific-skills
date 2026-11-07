# Pairwise



| Tool | Reference | Score | Best for |
|------|-----------|-------|----------|
| TM-align | Zhang & Skolnick 2005 NAR | TM-score | Single-chain pairwise; standard fold-similarity benchmark |
| US-align | Zhang et al 2022 Nat Methods | TM-score | Multi-chain protein, RNA, DNA, complexes; original successor to TM-align |
| Foldseek-Multimer | Kim et al 2025 Nat Methods | TM-score | Multi-chain complex search over thousands of complexes and up (the paper reports 3-4 orders of magnitude faster than US-align at database scale); no gain to expect for a single pair |
| CE | Shindyalov & Bourne 1998 Prot Eng | CE score | Combinatorial extension of fragments; PDB legacy; run as PyMOL `cealign` (Visualisation) |
| DALI | Holm 2022 NAR | Z-score | Distance-matrix alignment; superior at TM 0.3-0.5 (twilight fold); operates on AFDB at scale; local DaliLite in "Foldseek vs DALI" |
| Bio.PDB.Superimposer | Cock et al 2009 Bioinf | RMSD | Pure-Python superposition with known atom correspondence |

TM-align and US-align report two TM-scores by default: one normalised by chain 1 length (TM1) and one by chain 2 length (TM2); TM-align itself advises the reference (chain 2). Report both. For a whole-chain fold call use the smaller, `min(TM1, TM2)` (normalised by the longer chain). The larger only answers "is the smaller structure contained in the larger one" and gives false same-fold calls for short chains: with `max`, 3 of 86 different-fold pairs exceeded 0.5 (1PGA, 56 residues, against three unrelated kinases: 0.505-0.544); with `min`, none did (TM-align 20240303, 136 pairs of real PDB entries). The `-a T` flag adds a third score normalised by the average of the two chain lengths (useful for symmetric clustering); it needs the full output, because `-outfmt 2` rejects `-a`, `-u`, `-L` and `-d`. Thresholds are in the table below. RMSD alone is misleading: RMSD scales with length, depends on outliers, and an "optimization RMSD" (used to fit) is not the same as an "evaluation RMSD" (used to compare). Always report TM-score alongside RMSD and the aligned residue count (`(RMSD, n_atoms_used)`).

| Metric | Threshold | Source / Interpretation |
|--------|-----------|-------------------------|
| TM-score | > 0.5 | Same fold (Zhang & Skolnick 2004 Proteins; reported by TM-align, US-align, Foldseek, DALI) |
| TM-score | > 0.8 | Close topological match (not proof of homology) |
| TM-score | < 0.2 | Random structural similarity (not weak homology) |
| DALI Z-score | > 20 | Definitely homologous (Holm 2020 Methods Mol Biol) |
| DALI Z-score | 8 - 19 | Probable homology |
| DALI Z-score | 2 - 8 | Candidate; verify with TM-score or biology |
| DALI Z-score | < 2 | Not significant (random) |
| GDT-TS | > 50 | Correct fold (CASP/LGA scoring; reported by LGA, MaxCluster, OpenStructure and the Zhang-lab `TMscore` (below) -- not by TM-align/Foldseek) |
| LDDT | > 0.6 | Conventional cutoff for "correctly modelled" residue (Mariani et al 2013 Bioinf shows the score is fold-architecture-dependent and does not define a universal hard threshold; >= 0.6 is widely used in CAMEO and Foldseek output as a working cutoff) |
| RMSD | < 2 A over >100 residues | Strong superposition; below 1.5 A is excellent |

### Foldseek vs DALI: Modern Comparison

Both target structural homolog search, but with different strengths:

| Question | Foldseek | DALI |
|----------|----------|------|
| Find any same-fold homolog quickly | YES (3Di indexed; 1000-1M seq/s) | Slow (full distance matrix) |
| Sensitivity at TM > 0.5 (same-fold) | Comparable to TM-align | Slightly higher recall |
| Sensitivity at TM 0.3-0.5 (twilight fold) | Recall drops sharply | DALI Z-score retains signal (Holm 2022) |
| Alignment quality for indel-rich pairs | Local 3Di+AA; can miss large indels | Distance-matrix optimal handles indels well |
| AFDB-scale all-vs-all | Tractable | Now tractable post-2022 (Holm DALI server update) |

Practical workflow for remote homology: (1) Foldseek `easy-search` to retrieve same-fold candidates fast, (2) for hits with TM < 0.5 or with structural-indel rich alignments, re-align with DALI for higher-quality residue equivalences. The two tools are complementary; Foldseek-only misses some twilight-fold relationships, DALI-only misses the AFDB-scale opportunity.

Local DALI with DaliLite v5 (`import.pl`, `dali.pl`, the result-file trap and checked Z-scores) is in `references/dalilite.md`; the DALI web server is the alternative.

### TM-score Threshold Caveats

Apply the 0.5 same-fold rule only to globular domains of ~100-300 residues. Below ~60 residues do not call a fold from TM-score alone (a synthetic 10-residue helix scores 0.846 against myoglobin when normalised by its own length). Xu & Zhang 2010 give a length-aware significance for TM-score, but no tool here prints a p-value, so quote both TM-scores with the chain lengths.

### TM-align Pairwise Run

**Goal:** Compute TM-score and RMSD between two structures and write a superposed PDB.

**Approach:** Invoke TMalign or USalign with output flags; parse the structured output for downstream filtering.

```bash
TMalign chainA.pdb chainB.pdb -o sup            # writes sup.pdb (chain A superposed onto B, full atom) + sup*.pml
TMalign chainA.pdb chainB.pdb -outfmt 2

USalign chainA.pdb chainB.pdb -mol prot -outfmt 2
USalign complex_A.pdb complex_B.pdb -mm 1 -ter 0
```

`-outfmt 2` returns tabular output (one line per pair). `-o` takes a prefix and always appends `.pdb` (`-o sup.pdb` writes `sup.pdb.pdb`). For multi-chain complexes, US-align with `-mm 1 -ter 0` aligns full assemblies; TM-align is single-chain only. TM-align exits 0 with no data row for a chain of fewer than 3 residues or a file without protein atoms (`Sequence is too short <3!`, `Cannot parse file`), so check that a row came back (`examples/tm_align_pairwise.py` raises). US-align can segfault on a one-residue input.

**GDT-TS of a model against its native structure** (the `TMscore` binary from the same Zhang-lab source as TM-align):

```bash
TMscore model.pdb native.pdb -seq     # residues paired by sequence alignment; prints TM-score, RMSD, GDT-TS
```

Without `-seq`, TMscore pairs residues by residue number: an AlphaFold myoglobin model numbered one residue off from 1MBN gave GDT-TS 0.508 (TM 0.593) instead of 0.985 (TM 0.980, RMSD 0.705 over 153). Renumber, or use `-seq`.

**US-align multi-chain score interpretation.** US-align `-mm 1 -ter 0` prints two whole-complex TM-scores (normalised by structure 1 and by structure 2) and searches chain permutations itself. It also handles different stoichiometries directly: the hemoglobin dimer 1IRD against the tetramer 1A3N gives 0.977 (by the dimer) and 0.494 (by the tetramer), RMSD 0.94, 286 aligned residues. Report both: high by the smaller and about 0.5 by the larger means the small complex sits inside the large one. Do not use `-mm 4` for this: it is multiple-structure alignment (structure_2 is ignored) and segfaulted on two complexes. Per-chain TM-scores serve as a sanity check (low per-chain + high complex TM = topology match without local fold conservation, suggests promiscuous interaction not deep homology). See Zhang et al 2022 Nat Methods for the full mode taxonomy.

### Foldseek-Multimer for Database-Scale Complex Search

For multi-chain complex search over thousands of complexes or more, US-align is too slow and Foldseek-Multimer (`foldseek easy-multimersearch`, `easy-multimercluster`) is the default. Commands, the `_report` columns, the clustering trap and the speed evidence are in `references/foldseek-multimer.md`. Pairwise complex pair on a few hundred targets: stay with US-align `-mm 1 -ter 0`.

### Bio.PDB.Superimposer

**Goal:** Superpose a known atom correspondence and compute RMSD without running an external aligner.

**Approach:** Use when residue equivalence is already established (e.g. same sequence, different conformations). For unknown correspondence, prefer TM-align / US-align. Pair CA atoms by (chain, residue number, insertion code) over standard residues of the first model. Selecting every atom named `CA` also picks up Ca2+ ions, and pairing by list position then matches ions to residues: apo/holo calmodulin (1CFD/1CLL) printed 13.38 A that way instead of the true 10.83 A over 144 pairs, with no warning because both files had 148 "CA" atoms.

```bash
python scripts/biopython_superimposer.py reference.pdb mobile.pdb   # writes mobile_superposed.pdb
```

From Python, `from biopython_superimposer import superpose_ca` and call `superpose_ca(reference_structure, mobile_structure)`; it returns the `Superimposer`, the pair count and the number of pairs with different residue names.

The script refuses two cases: fewer than half the shorter chain paired, or more than 20% of paired residue names differing (offset numbering; `1MBN` vs `1A3N` paired by number would otherwise print 7.5 A over 141 pairs). Modified residues written as HETATM (selenomethionine `MSE`, phosphorylated residues; 1-4% of residues in the six real entries checked) are skipped by `residue.id[0] == ' '`; they are missing from the pairing, not mis-paired.



Structural alignments are read by viewing the superposed structures, not the sequence text:

```bash
# PyMOL headless: -cq runs without GUI/quietly; -d passes commands directly
pymol -cq -d "load reference.pdb; load mobile.pdb; super mobile, reference; ray 800,600; png fig.png"
```

PyMOL `super` performs cycle-fitting for distantly related structures (better than `align` when sequence identity is low; use `align` for closer matches); `cealign` runs CE; `tmalign` (PyMOL plugin) wraps TM-align.
