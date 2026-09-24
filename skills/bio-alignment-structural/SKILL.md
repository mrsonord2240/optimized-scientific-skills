---
name: bio-alignment-structural
description: Align protein structures using Foldseek 3Di, TM-align, US-align, DALI, or Foldmason for structural MSA. Score and superpose backbone coordinates when sequence identity is below the twilight zone or remote-homology detection is required. Use when sequence MSA fails (<25% identity), when the dark proteome is the target, when AlphaFoldDB / ESM Atlas search is needed, or when structural superposition is the goal.
tool_type: mixed
primary_tool: Foldseek
license: MIT
author: GPTomics
---

## Version Compatibility

Reference examples tested with: Foldseek 8+, TM-align 20220412+, US-align 20231222+, Foldmason 1+, BioPython 1.83+, pymol-open-source 3.0+. Every command below was run on Foldseek 10.941cd33, TM-align 20240303, US-align 20241108, Foldmason 4.dd3c235, DaliLite v5, PyMOL 3.1.0, Biopython 1.88.

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `foldseek version`, `foldmason version` (`--version` is rejected), `TMalign -v`, `USalign -v`
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed package and adapt the example to match the actual API rather than retrying.

# Structural Alignment

**Install.** `conda install -c bioconda foldseek usalign foldmason` (Python: `pip install biopython`).
Do not install bioconda `tmalign`: its binary is a 2018 Fortran build that dies with `libgfortran.so.3`. Build TM-align (and `TMscore`, used for GDT-TS below) from the US-align source (github.com/pylelab/USalign): `g++ -static -O3 -ffast-math -o TMalign TMalign.cpp` (same for `TMscore.cpp`), or use `USalign a.pdb b.pdb -mol prot`, which reproduced TM-align's TM-scores and RMSD to the printed digits. Foldseek databases: `foldseek databases` lists the names (`Alphafold/UniProt50`, `Alphafold/Swiss-Prot`, `PDB`, `CATH50`, ...); download with `foldseek databases <name> /path/to/db tmp/` (AFDB and PDB are tens of GB). DaliLite v5 (GPL, source build): ekhidna2.biocenter.helsinki.fi/dali/DaliLite.v5.tar.gz. MUSTANG 3.2.3 is a source build (Konagurthu et al. 2006).

**"Align two protein structures"** -> Compute backbone-aware superposition and a fold-similarity score (TM-score, RMSD, or LDDT).
- CLI pairwise: `TMalign A.pdb B.pdb`, `USalign A.pdb B.pdb`
- CLI search at scale: `foldseek easy-search query/ AFDB result.m8 tmp/`
- CLI structural MSA: `foldmason easy-msa structures/*.pdb out tmp/`
- Python pairwise: `Bio.PDB.Superimposer`, or `subprocess` wrapping `TMalign` / `USalign` (see `examples/tm_align_pairwise.py`)
- Scripted molecular-graphics superposition: PyMOL `super`/`cealign`

**"Find structural homologs of an AlphaFold model"** -> Search a structure database by 3Di-encoded structural alphabet (Foldseek) or by full TM-align rotation (DALI, US-align).

## When to Use Structural Alignment

| Sequence identity | Recommended approach |
|-------------------|---------------------|
| >= 40% | Sequence DP (Bio.Align, BLASTP) is sufficient |
| 25-40% | Sensitive sequence (MMseqs2, jackhmmer, profile-profile HHsearch) |
| 15-25% | Profile-profile (HHsearch) OR Foldseek if structures available |
| < 15% (dark proteome / twilight zone) | Foldseek (3Di), TM-align, US-align |

Sequence alignment below 15% identity is statistically indistinguishable from random pairings. The exact twilight-zone cutoff is length-dependent: Rost 1999 (Prot Eng) showed the curve drops to 25% at length 80, 20% at length 250 -- short alignments need higher identity for the same statistical signal, so a 15-25% rule of thumb is shorthand for "twilight zone for proteins of typical domain size (~150-300 residues)". If reasonable structural models exist (PDB, AlphaFoldDB, ESMFold), structural alignment is far more reliable in this regime.

### Twilight-Zone Threshold Exceptions

For families with strong functional or structural constraint (ribosomal proteins, class I aminoacyl-tRNA synthetases with HIGH/KMSKS motifs, histone fold, cytochrome c with CXXCH, or any long alignment with well-distributed conservation), sequence MSA may remain reliable down to ~12-15% identity. Verify with a profile-profile method (HHsearch) before committing to structural alignment.

## Pairwise Structural Alignment Tool Selection

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
python examples/biopython_superimposer.py reference.pdb mobile.pdb   # writes mobile_superposed.pdb
```

From Python, `from biopython_superimposer import superpose_ca` and call `superpose_ca(reference_structure, mobile_structure)`; it returns the `Superimposer`, the pair count and the number of pairs with different residue names.

The script refuses two cases: fewer than half the shorter chain paired, or more than 20% of paired residue names differing (offset numbering; `1MBN` vs `1A3N` paired by number would otherwise print 7.5 A over 141 pairs). Modified residues written as HETATM (selenomethionine `MSE`, phosphorylated residues; 1-4% of residues in the six real entries checked) are skipped by `residue.id[0] == ' '`; they are missing from the pairing, not mis-paired.

## Structural Search at Scale: Foldseek

Run Foldseek for AlphaFoldDB-scale structural search; it indexes a 20-letter 3Di alphabet for thousand- to million-fold speedup over TM-align at comparable same-fold sensitivity.

```bash
# Search query structures against AFDB (default: --alignment-type 2)
foldseek easy-search query.pdb afdb_database result.m8 tmp/

# Refine top hits with full TM-align rotation (slower but global TM-score); filter on alntmscore, not E-value
foldseek easy-search query.pdb afdb_database result.m8 tmp/ --alignment-type 1

# All-versus-all clustering at TM > 0.5
foldseek easy-cluster structures/*.pdb cluster_result tmp/ --tmscore-threshold 0.5

# Custom output columns
foldseek easy-search query.pdb afdb_database result.m8 tmp/ \
    --format-output query,target,evalue,alntmscore,qtmscore,ttmscore,lddt,bits
```

| Foldseek `--alignment-type` | Algorithm | Use when |
|-----------------------------|-----------|----------|
| 0 | 3Di Gotoh-Smith-Waterman (local) | Not recommended; 3Di alone, no amino-acid signal |
| 1 | TMalign (global) | Refine top hits with full TM-score; slow. E-values are meaningless (0.88-0.99 for true homologs, alnTM 0.93): an `evalue < 1e-3` filter returns 0 hits, so filter on `alntmscore` (and `lddt`) |
| 2 | 3Di+AA Gotoh-Smith-Waterman (local) | Default; best speed-sensitivity tradeoff |

`--max-seqs` (default 1000; `examples/foldseek_search.py` uses 200) silently caps the rows per query: a printed count equal to the cap is a cap, not the hit total.

### pLDDT-Filtering AlphaFold Structures Before Foldseek

Mask residues with pLDDT < 70 before Foldseek indexing or search; their backbone coordinates encode as random 3Di letters and contaminate hits below TM ~ 0.4. AFDB clusters from Barrio-Hernandez et al 2023 are pre-filtered; ESMFold predictions need the same step.

```bash
# Mask low-pLDDT residues during database creation (B-factor column stores pLDDT)
# The *.pdb glob requires shell expansion; for Python subprocess calls use glob.glob()
# to expand the file list before passing as argv.
foldseek createdb --mask-bfactor-threshold 70.0 *.pdb afdb_masked
```

## Structural Multiple Sequence Alignment

| Tool | Reference | Best for |
|------|-----------|----------|
| Foldmason | Gilchrist et al 2026 Science 391:485 | Billion-protein-scale structural MSA on AFDB; the structural counterpart to MAFFT |
| MUSTANG | Konagurthu et al 2006 Proteins | Pure structural MSA via residue-residue equivalences; `mustang-3.2.3 -i a.pdb b.pdb c.pdb -o out -F fasta` writes `out.afasta` + `out.pdb` (3 myoglobins: 3-row FASTA) |

Not covered, because no runnable path could be verified: T-Coffee Expresso and 3D-Coffee (T-Coffee 12.00.7: Expresso needs a BLAST server; 3D-Coffee retries retired RCSB/wwPDB endpoints and falls back to `proba_pair`), PROMALS3D and mTM-align (web servers), FATCAT, ChimeraX (licence-gated download) and the protein-language-model aligners TM-Vec, vcMSA, DEDAL and pLM-BLAST (GPU and model weights). Use Foldmason or MUSTANG.

### Foldmason easy-msa

`foldmason easy-msa structures/*.pdb result tmp/ --refine-iters 100 --refine-seed 42 --report-mode 1` writes `result_aa.fa`, `result_3di.fa`, `result.nw` and `result.html`. Always set `--refine-seed` (or `--refine-iters 0`): unseeded refinement is random. Rows are chains, not files. Per-column LDDT is under the JSON key `scores` (`--report-mode 2`), not `per_column_lddt`. Details and the extraction snippet: `references/foldmason.md`.

## AlphaFold Integration

The "predict-then-align" workflow has become standard for remote-homology problems:

1. Search sequence with MMseqs2 / jackhmmer to seed an MSA.
2. Predict structure with AlphaFold2 (ColabFold), ESMFold, or AlphaFold3 -- see `structural-biology/modern-structure-prediction` for prediction workflows.
3. Search the predicted structure against AFDB or PDB with Foldseek.
4. Use Foldmason to derive a structure-aware MSA from the hits.
5. Refine sequence MSA using the structure-derived column equivalences.

Search AFDB clusters (~2.3 M Foldseek-derived clusters from Barrio-Hernandez et al 2023) as the canonical remote-homology entry point; supplement with the ESM Atlas (~600 M ESMFold metagenomic structures) when AFDB recall is insufficient. For pLDDT semantics and curated AlphaFoldDB entry handling, see `structural-biology/alphafold-predictions`.

## Visualisation and Inspection

Structural alignments are read by viewing the superposed structures, not the sequence text:

```bash
# PyMOL headless: -cq runs without GUI/quietly; -d passes commands directly
pymol -cq -d "load reference.pdb; load mobile.pdb; super mobile, reference; ray 800,600; png fig.png"
```

PyMOL `super` performs cycle-fitting for distantly related structures (better than `align` when sequence identity is low; use `align` for closer matches); `cealign` runs CE; `tmalign` (PyMOL plugin) wraps TM-align.

## Decision Tree by Goal

| Goal | First-line tool |
|------|-----------------|
| Two structures, known correspondence | `Bio.PDB.Superimposer` |
| Two structures, unknown correspondence | TMalign or USalign |
| Multi-chain complex, pairwise | USalign with `-mm 1 -ter 0` |
| Multi-chain complex, database search | `foldseek easy-multimersearch` (Foldseek-Multimer; `references/foldseek-multimer.md`) |
| Twilight-fold homology (TM 0.3-0.5) | DALI via DaliLite `dali.pl` (Z-score ranks low-similarity hits better than Foldseek; `references/dalilite.md`) |
| Structural homolog search at AFDB scale (single chain) | `foldseek easy-search` |
| All-vs-all clustering of structures | `foldseek easy-cluster --tmscore-threshold 0.5` |
| Multiple structure alignment, < 100 chains | MUSTANG or Foldmason `easy-msa` (`references/foldmason.md`) |
| Multiple structure alignment, > 1000 chains | Foldmason `easy-msa` |
| Distant homology with no structures | Predict with ColabFold/ESMFold first, then Foldseek |

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| TM-align prints `Sequence is too short <3!` or `Cannot parse file` and no data row, exit 0 | Chain of fewer than 3 residues, or a ligand-only PDB | Filter to protein chains; test for a data row, never the exit code |
| Foldseek `Input <path> does not exist`, exit 1 | Wrong database path | Check the path; `foldseek databases` lists downloadable names (`PDB`, `Alphafold/Swiss-Prot`, ...) |
| Foldseek 0 rows, exit 0 | One-residue or unparsable query, or filters too tight | Check the query file; under `--alignment-type 1` filter on `alntmscore`, not E-value |
| Bio.PDB `Fixed and moving atom lists differ in size` | Unequal atom lists | Pair CA atoms by (chain, number, icode) with `residue.id[0] == ' '` (see Bio.PDB.Superimposer) |
| TM-score normalised by a fixed length | `-L N` was passed | Drop `-L`; the default normalises by each chain length |
| Foldmason `structuremsa died`, exit 1 | Only one input structure | Give at least two structures |

## Reference Files

| File | Read when |
|------|-----------|
| `references/foldseek-multimer.md` | Searching or clustering multi-chain complexes with Foldseek-Multimer (commands, `_report` columns, speed evidence, reporting convention) |
| `references/dalilite.md` | Running DALI locally with DaliLite v5 (`import.pl`, `dali.pl`, result-file trap) |
| `references/foldmason.md` | Building a structural MSA with Foldmason (seeding, output files, per-column LDDT JSON) |

## Related Skills

- alignment/multiple-alignment - Sequence MSA when identity > 25%; complementary for hybrid tools
- alignment/pairwise-alignment - Sequence pairwise; use first to filter before structural alignment
- alignment/msa-parsing - Parse and analyze structural MSA output for downstream metrics
- alignment/msa-statistics - Apply per-column conservation to structurally-derived MSAs
- alignment/alignment-io - Read and convert structural MSA files for downstream tools
- alignment/alignment-trimming - Trim structural MSAs the same way as sequence MSAs
- structural-biology/modern-structure-prediction - Predict structures (AlphaFold2/3, ESMFold) used as input
- structural-biology/alphafold-predictions - Curated AFDB / pLDDT handling
- structural-biology/structure-navigation - Map alignment columns to PDB residues
- phylogenetics/modern-tree-inference - Trees from structural MSAs (Foldmason output works directly)

## References

- van Kempen M et al. 2024. Fast and accurate protein structure search with Foldseek. Nat Biotech 42:243-246.
- Kim W, Mirdita M, Levy Karin E, Gilchrist CLM, Schweke H, Soding J, Levy E, Steinegger M. 2025. Rapid and sensitive protein complex alignment with Foldseek-Multimer. Nat Methods 22:469-472.
- Zhang Y, Skolnick J. 2005. TM-align: a protein structure alignment algorithm based on the TM-score. NAR 33:2302-2309.
- Zhang C, Shine M, Pyle AM, Zhang Y. 2022. US-align: universal structure alignments of proteins, nucleic acids, and macromolecular complexes. Nat Methods 19:1109-1115.
- Holm L. 2020. Using Dali for protein structure comparison. Methods Mol Biol 2112:29-42 (Z-score interpretation thresholds).
- Holm L. 2022. Dali server: structural unification of protein families. NAR 50:W210-W215.
- Gilchrist CLM et al. 2026. Foldmason: multiple protein structure alignment at scale with 3Di. Science 391(6784):485-488.
- Barrio-Hernandez I et al. 2023. Clustering predicted structures at the scale of the known protein universe. Nature 622:637-645.
- Xu J, Zhang Y. 2010. How significant is a protein structure similarity with TM-score = 0.5? Bioinf 26:889-895.
