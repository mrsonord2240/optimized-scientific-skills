# Intron Retention and Microexon Detection

## Intron Retention: Canonical vs Detained vs Co-Transcriptional Unspliced

Three biologically distinct states all called "IR" by generic tools:

1. **Canonical RI (cytoplasmic, NMD-substrate often)**: mature polyadenylated mRNA carries the intron; usually PTC-bearing and NMD-targeted, sometimes encoding an alternative protein.
2. **Detained intron (DI)** (Boutz 2015 *Genes Dev*): nuclear-localized, mature transcripts retaining a specific intron; a regulated reservoir released into translation upon signaling.
3. **Co-transcriptional unspliced**: nascent pre-mRNA captured before splicing complete; not a regulated state.

**Library prep determines which state(s) are visible:**
- Poly(A) selection: enriches (1), depletes (2)/(3)
- rRNA depletion (cytoplasmic): captures (1)
- rRNA depletion (whole cell or nuclear): captures all three

**To distinguish DI from canonical RI:** subcellular fractionation (nuclear vs cytoplasmic RNA-seq), or NMD inhibitor (cycloheximide, NMDi-14) treatment — canonical RI mRNA increases under NMD inhibition; DI does not.

```bash
# One-time reference from a STAR genome index (IRFinder 1.3.1; BuildRef alone takes an Ensembl FTP URL)
IRFinder -m BuildRefFromSTARRef -r REF -x star_index/
# Quantify; single-end takes one FASTQ, paired-end two. Results: ir_output/IRFinder-IR-nondir.txt
IRFinder -m FastQ -r REF -t 8 -d ir_output R1.fastq.gz [R2.fastq.gz]
```

Executed with IRFinder 1.3.1. IRFinder-S (Lorenzi 2021 *Genome Biol*) adds CNN-based filtering of true IR vs noise and was not executed here. iREAD and S-IRFindeR (Broseus & Ritchie 2020 *bioRxiv*) are alternatives, also not executed.

## Microexon Detection

Microexons (3-27 nt, neural-enriched, SRRM4-regulated; Irimia 2014 *Cell*) are missed by default short-read aligners requiring 20-30 nt anchors. Options:

| Approach | Tool | Notes |
|----------|------|-------|
| Curated database lookup | VAST-TOOLS + VastDB | Cross-species, microexon-aware (Tapial 2017 *Genome Res*) |
| De novo discovery | MicroExonator (Parada 2021 *Genome Biol*) | Snakemake pipeline |
| Tune the upstream aligner | `STAR --alignSJoverhangMin 6 --alignSJDBoverhangMin 1 --outFilterMismatchNoverReadLmax 0.04` | rMATS itself cannot recover microexons that STAR didn't pass through; lower DB-junction overhang to 1 (trusts annotated microexon coords) and combine with strict mismatch filter. Typical AS pipelines use STAR 8/3 which is too strict for microexons |
| Long-read sequencing | PacBio Iso-Seq, ONT | Solves the problem entirely; reads span microexons fully |

For brain / neural tissue or autism-spectrum studies, **microexon analysis must be explicit** — default short-read pipelines underdetect them (~70% reported; not reproduced here).

