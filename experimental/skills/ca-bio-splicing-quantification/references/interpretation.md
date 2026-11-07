## Quality Thresholds

| Metric | Threshold | Source / Rationale |
|--------|-----------|---------------------|
| Junction reads per replicate | >=10-20 (per-replicate minimum) | Empirical PSI variance becomes <0.05 above this; below, PSI becomes a coin flip |
| PSI dynamic range | mean PSI 0.05-0.95 | Outside is near-constitutive; `filter_reliable_events` applies this (`psi_range`); not a rMATS or SUPPA2 default |
| Missing values | <50% of samples | Higher missingness indicates low expression — re-test with subset |
| Read length | >=75nt PE preferred; >=100nt for microexons | Shorter single-end reads bias junction detection toward shorter exons (convention) |
| Library | rRNA depletion for IR analysis; poly(A) acceptable for cassette | poly(A) selection loses pre-mRNA/intronic signal (convention) |
| STAR 2-pass | Cohort-style preferred over per-sample basic | Veeneman 2016 *Bioinformatics*: >=94% novel junction recovery |
| MAJIQ minreads / minpos | --minreads 10 --minpos 3 | Default; lower for low-coverage |
| leafcutter -m | 50 reads per cluster | Higher for rare events; lower for sQTL discovery |
| Anchor length | >=8 nt for short-read | Below this, false-positive junctions dominate (CIGAR-N noise) |

## Decision Tree by Scenario

| Scenario | Recommended tool(s) | Why |
|----------|----------------------|-----|
| Standard cassette analysis, n>=3, GENCODE-annotated | rMATS-turbo + leafcutter (concordance) | Default workflow; complementary algorithmic families |
| Non-model organism, no GENCODE-grade annotation | leafcutter + de novo discovery | Annotation-free |
| Heterogeneous cohort, n>=10 vs n>=10 (clinical, GTEx-style) | MAJIQ V3 with HET module | HET designed for between-sample variability dominance |
| Low coverage / few replicates (n=2 vs n=2) | Shiba | Junction-imbalance correction; SOTA at low coverage in 2025 benchmarks |
| Cross-species comparative (vertebrate panel) | VAST-TOOLS + VastDB | Orthology-aware events; ExOrthist co-tool |
| TPM-only available (no BAMs) | SUPPA2 | Annotation-bound but fast |
| Microexon focus (neural / ASD) | VAST-TOOLS or MicroExonator | Default tools systematically miss microexons |
| Intron retention focus | IRFinder (rRNA-depleted library) | Coverage-aware |
| Detained introns specifically | IRFinder + nuclear/cytoplasmic fractionation | Required to separate DI from cytoplasmic RI |
| Long reads available | rMATS-long, FLAIR, IsoQuant | Full-isoform resolution; see bio-long-read-splicing |
| Single-cell (full-length plate) | MARVEL, BRIE2 | See bio-single-cell-splicing |
| Single-cell (10X 3') | Likely don't attempt; consider Sierra for APA | 10X 3' chemistry insufficient for AS |

## Output Interpretation

PSI ranges 0 to 1: 1 = always included, 0 = always skipped, 0.5 = balanced. Sign of `IncLevelDifference` matches `--b1 minus --b2` group order — always document which is which in publications.

**NMD direction matters:** an increase in PSI of a poison exon (PTC-introducing) **decreases** functional protein due to NMD. Always check whether the alternative form is PTC-bearing using ORF-aware annotation (IsoformSwitchAnalyzeR consequences, or manual stop-codon distance check vs last exon-exon junction).

**Disease signatures:**
- **SF3B1** mutations (MDS, CLL, uveal melanoma): cryptic 3'ss ~10-30 nt upstream of canonical (Darman 2015 *Cell Rep*). Look for clustered A3SS hits.
- **U2AF1** mutations (lung adeno, MDS): altered preferences at 3'ss -3 position; cassette-exon shifts.
- **TDP-43 loss** (ALS/FTD): de novo cryptic exons in UNC13A, STMN2, ATG4B (Brown 2022 *Nature*; Klim 2019 *Nat Neurosci*) — annotation-free tools required (leafcutter denovo).

## Common Pitfalls

- Treating AFE/ALE as splicing — these are typically promoter-driven (AFE) or APA-driven (ALE), not spliceosomal. Confirm with FANTOM CAGE or 3'-end-seq.
- Confusing detained introns with NMD-targeted RI — both call as "IR" but have opposite biological fates.
- Using poly(A) libraries for IR analysis — biases toward mature transcripts, depletes pre-mRNA.
- Single-end short reads — junction-spanning reads need >=8nt overhang on both sides; biases toward shorter exons.
- Quoting "PSI of the gene" from MAJIQ LSV output — only per-junction PSI within an LSV is meaningful.
- Skipping STAR 2-pass — loses novel junctions (~14% reported; not reproduced here); matters for any non-canonical organism or condition.
- Trusting rMATS calls without `--novelSS` when STAR 2-pass found new junctions — rMATS will only quantify pre-annotated events.
