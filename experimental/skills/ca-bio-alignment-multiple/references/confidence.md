# Confidence



**Goal:** Quantify alignment uncertainty by generating multiple HMM-perturbed alignments and measuring column consistency.

**Approach:** MUSCLE5 (Edgar 2022 Nat Comm) ships two ensemble modes: `-stratified` (16 replicates by default: the `-replicates` flag defaults to 4 HMM-perturbation seeds x 4 guide-tree permutations) and `-diversified` (100 replicates by default). Both belong to the `-align` (PPP) command, not `-super5` -- `-super5` only accepts `-perm`/`-perturb` for a single replicate and rejects `-stratified`/`-diversified`/`-replicates` outright (checked on MUSCLE 5.3.win64: each fails with "`<flag>` not supported"). Both write an Ensemble FASTA (.efa) file containing all replicates; column-level confidence is the fraction of replicates that place a given residue pair in the same column. `-perturb SEED` is a separate flag that sets the HMM-perturbation random seed, not an ensemble selector.

```bash
# Stratified ensemble: 16 replicates (4 HMM-perturbation seeds x 4 guide-tree permutations)
# Verified on MUSCLE 5.3.win64 / prot15_unaligned.fa (15 seqs): produced 240 aligned records = 16 x 15
muscle -align input.fasta -stratified -output ensemble.efa

# Diversified ensemble: 100 replicates exploring guide-tree and HMM space
# Verified: produced 1500 aligned records = 100 x 15
muscle -align input.fasta -diversified -output ensemble.efa

# Optional: change replicate count and HMM-perturbation seed
muscle -align input.fasta -stratified -replicates 8 -perturb 42 -output ensemble.efa
```

`-align` scales to peak accuracy only up to ~1000 sequences. Above that you need `-super5`, which cannot write `.efa` directly. Run it once per replicate with a distinct `-perm`/`-perturb` pair, then combine the single alignments into one ensemble with `-fa2efa` (verified: 3 `-super5` single-replicate runs combined into one 45-record `.efa` = 3 x 15 seqs):

```bash
# >1000 sequences: build replicates individually with -super5, then merge into one ensemble
muscle -super5 input.fasta -perm none -perturb 1 -output rep1.afa
muscle -super5 input.fasta -perm abc  -perturb 2 -output rep2.afa
muscle -super5 input.fasta -perm acb  -perturb 3 -output rep3.afa
printf 'rep1.afa\nrep2.afa\nrep3.afa\n' > replicates.txt
muscle -fa2efa replicates.txt -output ensemble.efa
```

The .efa output is consumed downstream to derive confidence-weighted bootstrap support: each replicate is fed to a tree builder and the resulting trees combined (Edgar 2022 supplement). Columns consistently aligned across replicates are reliable; high-divergence regions diverge between replicates and should be flagged before phylogenetic inference.



For publication-grade phylogenetics or selection analysis, alignment uncertainty must be quantified per column and unreliable columns excluded BEFORE downstream inference -- not after.

| Method | Tool | Output |
|--------|------|--------|
| Guide-tree perturbation | GUIDANCE2 (Sela et al 2015 NAR) | Per-column and per-residue reliability score (0-1); default cutoff 0.93. No GUIDANCE3 has been released; deep-learning successors are exploratory only. |
| Library-consistency | T-Coffee TCS (Chang et al 2014 MBE) | Per-column reliability score via `-evaluate -output score_ascii` |
| HMM ensemble | MUSCLE5 `-stratified` / `-diversified` | EFA file; column confidence = fraction of replicates supporting it |
| Co-optimal alignments | HoT (Landan & Graur 2007 MBE) | Heads-or-tails alternate optimal alignment for each column |

Mask columns below the reliability threshold before phylogenetic inference. Trees built from filtered alignments are markedly more stable to method choice.

**GUIDANCE2 thresholds are tool-specific.** The 0.93 default cutoff is calibrated against MAFFT-LINSI in Sela et al 2015. Running GUIDANCE2 on top of MUSCLE5 or PRANK uses different perturbation distributions and produces scores on a slightly different scale. For non-MAFFT base aligners, calibrate the threshold empirically using a benchmark with known-correct columns or use the tool-native ensemble scoring (MUSCLE5 `-stratified`) instead.

## Post-Alignment Validation Checklist

Before proceeding to downstream analysis, verify alignment quality:

1. **Visual inspection**: Scan for columns of mostly gaps with scattered residues (hallmark of misalignment)
2. **Gap distribution**: A high gap fraction (>50% of columns with gaps) is a prompt to inspect further, not a verdict on its own -- a correct indel-rich alignment (true divergent-family history) can clear this threshold too (audited case: 55.7% gapped columns on a genuinely correct alignment, SP 0.966 vs truth). Read it alongside the homology/orientation screen below, not instead of it.
3. **Sequence identity**: If average pairwise identity is <25% for proteins, alignment reliability is questionable
4. **Outlier sequences**: Sequences with excessive gaps relative to others are a lead to check, not a homology verdict -- a contaminant's gap fraction can sit inside the range of the real homologs (audited case: a non-homologous contig at 0.18 next to homologs at 0.12-0.16, indistinguishable by gap fraction alone) while scoring far below them on a direct homology/orientation screen. Run the screen in "When NOT to Run MSA" below rather than deciding from gap counts alone.
5. **Conservation pattern**: Functional domains should show clear conservation; absence of expected conserved motifs suggests alignment error or non-homology
6. **Run GUIDANCE2 or MUSCLE5 ensemble**: Quantify alignment confidence per column before phylogenetic inference
