# Trimming, Gap Treatment and Unreliable Regions

Read when the request is about trimming tools, how gaps enter a phylogenetic analysis, or masking unreliable columns (MUSCLE5 column confidence).

## Alignment Trimming

Trimming controversy and tool selection (ClipKIT, trimAl, BMGE, Divvier, HMMcleaner, Noisy) is the subject of a dedicated skill. Use this short decision matrix for routing:

| Goal | First-line tool |
|------|-----------------|
| Phylogenetic-tree input | ClipKIT `kpic-smart-gap` (Steenwyk et al 2020 PLOS Bio) |
| HMM profile building | trimAl `-gappyout` (Capella-Gutierrez et al 2009 Bioinf) |
| Selection / dN/dS input | Avoid aggressive trimming; mask by MUSCLE5 ensemble column confidence (see Identifying Unreliable Alignment Regions) |
| Deep prokaryotic phylogenomics | BMGE (Criscuolo & Gribaldo 2010 BMC Evol Biol) |
| Preserve column-mapping for residue-level analysis | trimAl `-colnumbering` |

Prefer ClipKIT `kpic-smart-gap` over traditional gap-only removal for trees, and note that aggressive trimming (>20-30% of sites) can hurt tree quality. See alignment/alignment-trimming for full mode comparisons, decision trees, and runnable examples.

## Gap Handling for Phylogenetics

How gaps are treated in downstream phylogenetic analysis significantly affects tree topology:

| Treatment | Method | Tradeoff |
|-----------|--------|----------|
| Missing data (default) | Gaps = unknown character | Most common; can be statistically inconsistent under ML |
| Fifth state | Gap = 5th nucleotide | Biologically problematic (gaps of different lengths treated equally) |
| Simple indel coding | Each unique indel coded as binary character | Most biologically realistic; adds phylogenetic signal |

For slow- to mid-rate datasets where indels are phylogenetically informative, prefer SIC indel coding (use fifth-state only as a sensitivity check); for rapidly-evolving datasets (intra-species, ITS regions, retroelement-rich plant genomes), default to missing-data treatment because gap homology is unreliable. Run a sensitivity analysis comparing treatments before drawing topological conclusions.

## Identifying Unreliable Alignment Regions

Columns exhibiting **both** high gap fraction AND low conservation are the strongest indicators of alignment uncertainty. These often reflect guide tree artifacts rather than true evolutionary events; divergent sequences disproportionately introduce gaps. Before phylogenetic analysis:

1. Flag columns with gap fraction >50%, which may be alignment artifacts
2. Check if gappy regions coincide with insertions in a single divergent sequence (remove that sequence and re-align)
3. For critical analyses, get per-column confidence from a MUSCLE5 ensemble and mask low-confidence columns (checked on MUSCLE 5.3; needs the unaligned sequences, not the alignment). GUIDANCE2 is not offered here: its stand-alone package is no longer downloadable.

```bash
muscle -align seqs.fa -stratified -output ens.efa      # ensemble of replicate alignments
muscle -maxcc ens.efa -output maxcc.afa                # stderr line "CC min .., best <name>", e.g. acb.2
muscle -addconfseq ens.efa -output ens_cc.efa          # adds per-column confidence (CC) rows
python examples/muscle5_column_confidence.py ens_cc.efa acb.2 0.9 masked.fa
```

Run the script without the last two arguments first: it prints how many columns survive at CC 0.5/0.7/0.9/0.99, and there is no calibrated cut-off (GUIDANCE2's 0.93 does not transfer). On 8 UniProt globins 11 of 155 columns had CC < 0.9, including every column whose residue pairing differed in more than 10% of the 16 replicates (independent replicate-agreement check).
