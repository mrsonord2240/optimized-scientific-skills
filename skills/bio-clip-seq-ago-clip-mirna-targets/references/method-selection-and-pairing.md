# Method selection and pairing rules

## Methods taxonomy

| Method | Year | miRNA-target pairing | Chimera enrichment | Strength | Fails when |
|---|---:|---|---|---|---|
| HITS-CLIP for AGO | 2009 | Indirect, computational seed assignment | None | Original and widely cited | A direct miRNA assignment is required |
| PAR-CLIP for AGO | 2010 | Indirect | None | T-to-C signature at the crosslink position | Cells do not tolerate or incorporate 4SU |
| CLASH | 2013 | Direct chimera | None | Pan-Argonaute interactome | Chimera recovery is too low |
| AGO-CLEAR-CLIP | 2015 | Direct chimera | None, incidental recovery | Direct miRNA-target pairs | Library depth is inadequate for a low chimera rate |
| HEAP | 2020 | Indirect, with a chimeric step | None | HaloTag-Ago2 profiling in vivo | The system is not the transgenic mouse model |
| Chimeric eCLIP / miR-eCLIP | 2022 | Direct chimera | Probe or PCR enrichment | Deep profiling, including selected miRNAs | Specialized library preparation is unavailable |
| AGO-IP-microarray | 2007 | Indirect | None | Historical predecessor to AGO CLIP | Crosslink-dependent or transient-site resolution is required |

Methodology evolves. Verify current primary literature and the exact library
protocol before describing one method as canonical or state of the art.

## Direct versus computational pairing

Chimeric methods covalently join the miRNA and target fragment during library
preparation. The two aligned portions therefore identify the pair directly,
including interactions that lack a canonical seed. Standard AGO CLIP enriches
AGO-bound RNA but does not preserve which miRNA occupied each site. Seed scans
then infer candidate miRNAs and may assign many candidates to one peak.

The inherited source material includes numeric planning heuristics for
unenriched chimera recovery and sequencing depth without claim-level
provenance. Do not use those numbers as quality thresholds. Base depth and
enrichment decisions on the exact protocol's primary source, controls, and
protocol-matched pilot recovery.

## Canonical seed classes

| Seed type | Source pairing definition | Position 1 | Interpretation |
|---|---|---|---|
| 8mer | miRNA positions 2-8 plus target A opposite position 1 | A required | Strong canonical site |
| 7mer-m8 | miRNA positions 2-8 | Any | Strong canonical site |
| 7mer-A1 | miRNA positions 2-7 plus target A opposite position 1 | A required | Canonical site |
| 6mer | miRNA positions 2-7 | Any | Weak and common |
| 6mer-A1 | miRNA positions 2-6 plus target A opposite position 1 | A required | Weak |
| 3'-compensatory | Weak seed plus strong pairing around miRNA positions 12-17 | Any | Can be missed by seed-only assignment |
| Central pairing | Central pairing without a seed | Any | Rare; investigate separately |

For TargetScan integration, use a versioned prediction set, retain 7mer and
8mer classes for the primary high-confidence overlap, and report 6mer sites
separately if included.

## Decision guide

| Scenario | Preferred source route | Reason |
|---|---|---|
| Modern direct target discovery | Chimeric eCLIP / miR-eCLIP | Direct chimera evidence and optional enrichment |
| Deep target list for a specific miRNA | miR-eCLIP with a matching probe | Enrichment improves per-miRNA coverage |
| Novel or non-canonical interactions | CLEAR-CLIP or chimeric eCLIP | Does not require a canonical seed prior |
| Initial AGO-binding discovery | AGO HITS-CLIP or eCLIP | Cost-effective site mapping |
| Cross-species comparison | TargetScan plus AGO CLIP in each species | Conservation prediction plus species-matched binding |
| miRNA perturbation | Knockdown/knockout plus AGO CLIP and differential analysis | Measures perturbation-associated changes |
| Single-prediction validation | miR-eCLIP for the miRNA plus an orthogonal reporter assay | Direct and functional evidence |

## Evidence reconciliation

| Pattern | Likely explanation | Interpretation or action |
|---|---|---|
| Chimera present, TargetScan site absent | Non-canonical or 3'-supplementary pairing | Retain the direct observation; characterize the duplex |
| TargetScan site and AGO peak, no chimera | Functional site with stochastic or shallow chimera capture | Report as AGO-supported prediction |
| TargetScan site, no AGO peak | Assay/context may lack detectable AGO occupancy | Report as unsupported in that assay and context, not disproven |
| AGO peak, no TargetScan match | Non-canonical site or unmodeled expressed miRNA | Investigate without forcing an assignment |
| Per-miRNA counts differ by orders of magnitude | Expression and library recovery differ | Normalize interpretation to matched expression and depth |
| Hyb chimera fraction near 1% | May reflect protocol, depth, mapping, or recovery | Compare to matched controls and the exact source protocol before redesign |
| Chimera tools give different counts | Algorithm sensitivity and ambiguity policies differ | Compare versioned tools and report discordant assignments explicitly |
| HEAP and eCLIP disagree | Species, tissue, and in-vivo/cell-line context differ | Preserve both contexts rather than merging them |

## Source references

- Chi SW et al. 2009. *Nature* 460:479. AGO HITS-CLIP.
- Hafner M et al. 2010. *Cell* 141:129. PAR-CLIP for AGO.
- Helwak A et al. 2013. *Cell* 153:654. CLASH.
- Travis AJ et al. 2014. *Methods* 65:263. Hyb pipeline.
- Moore MJ et al. 2015. *Nature Communications* 6:8864. CLEAR-CLIP.
- Li X et al. 2020. *Molecular Cell* 79:167. HEAP.
- Agarwal V et al. 2015. *eLife* 4:e05005. TargetScan 7.0.
- Lewis BP et al. 2003. *Cell* 115:787. Canonical seed rules.
- Bartel DP. 2018. *Cell* 173:20. miRNA targeting principles.
- McGeary SE et al. 2019. *Science* 366:eaav1741. Quantitative target prediction associated with TargetScan 8.0.

## Related source skills

- `clip-seq/clip-peak-calling`
- `clip-seq/binding-site-annotation`
- `clip-seq/clip-motif-analysis`
- `clip-seq/differential-clip`
- `clip-seq/m6a-clip`
- `small-rna-seq/target-prediction`
- `small-rna-seq/differential-mirna`
- `small-rna-seq/mirdeep2-analysis`
