# Specialist Trimmers: HMMcleaner, Gblocks, PhyIN

Moved out of `SKILL.md` (bio-alignment-trimming) because these three are needed less often than
ClipKIT, trimAl, BMGE and Divvier -- consult this file when the Goal-Driven Tool Selection table or
Decision Tree in `SKILL.md` points you here.

## HMMcleaner: Per-Residue Cleaning

Run HMMcleaner before column trimming when cross-contamination or annotation errors are suspected; it masks suspect residues with `X` rather than removing columns. Apply only to alignments with >= 15 sequences (below that the per-residue pHMM has too little signal and over-flags divergent residues).

```bash
HmmCleaner.pl input.fasta
```

Output is `input_hmm.fasta` with low-confidence residues masked. Used in OrthoMaM and Compositional Heterogeneity-aware phylogenomic pipelines where cross-contamination from sequencing or annotation errors is suspected.

**Sample-size sensitivity:** HMMcleaner builds a per-residue HMM from the OTHER sequences at each position, so on small alignments (<15 sequences) the HMM has insufficient signal and over-flags real divergent residues as contamination. Di Franco et al 2019 show two related effects: pHMMs built from few sequences carry less signal, lowering HMMcleaner sensitivity, and its false positives are driven mainly by evolutionary divergence -- specificity falls with gap frequency, evolutionary rate, and the fraction of ambiguously-aligned regions rather than by a fixed sequence count (global specificity ~94% at default settings). OrthoMaM v11 and PhyloHerb pipelines apply HMMcleaner only to alignments with >= 15 sequences. For smaller alignments, manual inspection or BLAST-based contamination screening is more reliable.

## Gblocks: The Legacy Default

Use Gblocks only when matching a legacy pipeline; prefer ClipKIT or trimAl for new work. Default parameters are too aggressive -- relax them. `-b1` and `-b2` are integer sequence COUNTS, not percentages: `-b1` must be greater than half the number of sequences N (default N/2 + 1), and `-b2` at least `-b1` (default 85% of N). Relaxed example for N = 20 sequences:

```bash
Gblocks input.fasta -t=p -b1=11 -b2=11 -b3=10 -b4=5 -b5=h
```

## PhyIN: Phylogenetic Incompatibility Trimming

PhyIN (Maddison 2024 PeerJ) is a complementary trimmer that flags neighbouring columns whose split patterns are phylogenetically incompatible (i.e. cannot share a tree). It targets a different failure mode than gap-based or entropy-based trimmers: ClipKIT, trimAl, and BMGE all evaluate columns independently, so a region with consistent gap structure passes their filters even if its character patterns are pairwise tree-incompatible (alignment artefact). PhyIN catches that case.

```bash
# PhyIN v1.0 is a single Python script; DNA/RNA FASTA alignments only
python phyin.py -input input.fasta -output trimmed.fasta -b 10 -d 2 -p 0.5
```

`-b` block length over which conflict is assessed, `-d` neighbour distance surveyed, `-p` proportion of conflicting neighbours that triggers removal; `-e` / `-not.e` treat gaps as an extra state or not. Use PhyIN as a SECOND-PASS trimmer after ClipKIT/trimAl when alignment artefact is the suspected source of incongruence in a phylogenomic dataset; not a replacement for first-pass column filtering. PhyIN's incompatibility test is signal-direction agnostic, so it does not distinguish "the alignment is wrong here" from "true incongruent locus" (e.g. introgression, ILS); only gene-tree comparison after tree-building can disambiguate. Reference: PeerJ 2024 paper.

## References

- Maddison WP. 2024. PhyIN: trimming alignments by phylogenetic incompatibilities among neighbouring sites. PeerJ 12:e18504.
