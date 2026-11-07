# Clipkit



Pick a `-m` mode by what the downstream task needs the alignment to KEEP. Run `clipkit --help` for the full list of 15 modes.

| Pick this mode | When |
|----------------|------|
| `smart-gap` | Default for single genes and supermatrices (ClipKIT's own default); keeps singleton and constant columns, so branch lengths are preserved |
| `kpic-smart-gap` | Topology-focused option for taxonomically balanced datasets only. Drops all singleton columns (typically 40-55% of a gene). Always compare outgroup/stem branch lengths before vs after; do not use when branch lengths feed dating, rate estimates or rooting |
| `gappy` / `gappyout` | Need an explicit gap threshold or trimAl-like behaviour |
| `kpi-smart-gap` / `kpi` | Parsimony-only input. Drops constant AND singleton columns, which inflates ML branch lengths; do not use for ML/Bayesian trees |
| `c3` / `cst` / `entropy` / `block-gappy` / `composition-bias` / `heterotachy` | Specialist needs: third-codon stripping, manual masks, entropy thresholding, conserved-block preservation, long-branch mitigation, or heterotachy filtering -- consult `clipkit --help` |

**Goal:** Run ClipKIT on an MSA with the mode appropriate to the downstream goal.

**Approach:** Pick a `-m` mode from the table above, optionally enable `--log` for reproducibility, and write the trimmed alignment in the desired format.

```bash
clipkit input.fasta -m smart-gap -o trimmed.fasta
clipkit input.fasta -m gappy -g 0.9 -o trimmed.fasta
clipkit input.fasta -m smart-gap -of phylip -o trimmed.phy
clipkit input.fasta -m smart-gap --log -o trimmed.fasta
```

The `--log` flag writes `<output>.log` with a per-column kept/removed record -- essential for reproducibility audits. ClipKIT also has a Python API (`clipkit.api.clipkit`) for in-memory trimming.

### ClipKIT kpic Failure Modes

`kpic` keeps parsimony-informative sites ("at least 2 sequences with each of at least 2 different residues") and constant sites, and removes singletons. On taxonomically unbalanced datasets (e.g. many close relatives + a few outgroups), a column where a single outgroup carries a private state is a singleton and is removed, so much of the outgroup-vs-ingroup signal is lost. In a simulated 30-strain + 3-outgroup DNA alignment, `kpic-smart-gap` removed most such columns: topology was unchanged, but the ingroup/outgroup stem length went from 0.44 to 1.2 substitutions/site (true 0.44) and one outgroup's terminal branch collapsed to 0. `smart-gap` left every branch length unchanged. Mitigations:
- Use `smart-gap` on unbalanced datasets and whenever branch lengths matter (dating, rates, rooting)
- If `kpic-smart-gap` is used, compare outgroup stem and terminal branch lengths before vs after trimming; if they move, discard the kpic result
- ClipKIT's site classification is count-based and unweighted (no sequence weighting by clade); manual subsampling to balance clades is the practical fix

BMGE's `-w` option is a sliding window for smoothing entropy values (default 3), not a clade-aware weighting.
