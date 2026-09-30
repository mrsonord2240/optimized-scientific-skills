---
name: bio-atac-seq-nucleosome-positioning
description: Map nucleosome center positions, occupancy, and fuzziness from ATAC-seq fragment-size patterns using NucleoATAC, ATACseqQC, DANPOS3, or scprinter. Use when characterizing nucleosome organization at promoters and enhancers, calling +1/-1 nucleosomes flanking NFRs, generating V-plots for chromatin structure visualization, or comparing nucleosome positioning between conditions.
license: MIT
tool_type: mixed
primary_tool: NucleoATAC
category: Data Analysis
author: GPTomics
---

# Nucleosome positioning from ATAC-seq

Tn5 cuts twice in naked accessible DNA (short fragments) and once on each side of a nucleosome (fragments of about nucleosome plus linker). Use fragment-size classes to call nucleosome centers, occupancy, and fuzziness, and to read spacing around regulatory elements.

## Inputs

Deduplicated, MAPQ-filtered, chrM-stripped paired-end BAM with about 30M or more nuclear reads. Mono-nucleosome analysis requires paired-end data. Use a genome build matched to the TSS or feature BED.

## Workflow

1. Confirm the BAM is paired-end and the feature coordinates (TSS, motif centers, peaks) match its genome build.
2. Plot fragment-size distribution and a V-plot at TSSs first ([`scripts/vplot.py`](scripts/vplot.py); R alternative in [`scripts/nucleosome_analysis.R`](scripts/nucleosome_analysis.R)). Give `vplot.py` a strand-column (6) BED of TSSs so minus-strand features are mirrored. A classic V or W means positioning is recoverable; a flat horizontal band means it is not, so stop and report the diagnosis instead of calling nucleosomes.
3. Choose a tool by goal:
   - per-base occupancy and nucleosome calls: NucleoATAC (`nucleoatac run --bed regions.bed --bam sample.bam --fasta genome.fa --out out/ --cores N`) or scprinter (no tested workflow in this Skill; follow its own documentation);
   - V-plot and fragment-class views in R: ATACseqQC;
   - differential positioning between conditions: DANPOS3 `dpos`;
   - single-cell: scprinter;
   - NRL estimate: [`scripts/estimate_nrl.py`](scripts/estimate_nrl.py).
4. Use regions of 500 bp or more (merge adjacent peaks with bedtools) for NucleoATAC.
5. For +1 nucleosomes, slop TSS intervals (`bedtools slop -l 200 -r 1000`), run NucleoATAC, and take the first called nucleosome downstream of each TSS from the nucpos BED. Restrict to protein-coding TSSs in the matched build.
6. Use the fragment windows NFR under 100 bp, mono 180-247 bp, di 315-473 bp. Do not filter mono-nucleosomes at 147 bp; the fragment includes linker. Adjust for the organism's NRL.
7. For DANPOS3 differential calls, require `abs(treat2control_dis)` of at least 30 bp and `point_diff_FDR` under 0.05 (filter command in the method reference).
8. Record tool versions, region definition, fragment windows, and thresholds.

## Routed material

- [`references/method-reference.md`](references/method-reference.md): NRL by organism, fragment classes, V-plot patterns, tool taxonomy and failure modes, DANPOS3 recipes, fuzziness bands, histone-variant and long-read notes, common errors, citations. Read it when interpreting patterns, choosing or troubleshooting a tool, or setting DANPOS3 parameters.
- [`references/usage-guide.md`](references/usage-guide.md): installation commands and example requests.
- Scripts: [`scripts/vplot.py`](scripts/vplot.py) (`vplot.py sample.bam features.bed [out.png]`), [`scripts/estimate_nrl.py`](scripts/estimate_nrl.py) (`estimate_nrl.py sample.bam`; reports the 150-250 bp mono-fragment mode or exits with a no-periodicity message), [`scripts/nucleosome_analysis.R`](scripts/nucleosome_analysis.R) (`Rscript nucleosome_analysis.R sample.bam [prefix] [tss.bed]`; hg38 ATACseqQC pipeline).

## Caveats

- NucleoATAC 0.3.4 is unmaintained and runs only on Python 2.7; see the method reference for its install route and the required nucpos check.
- Pioneer-factor asymmetry in V-plots is biological signal, not an artefact.
- ATAC alone cannot call H2A.Z or other variant nucleosomes; fragment size only generates hypotheses.

## License and provenance

This derived Skill retains GPTomics/bioSkills material by Domen Jemec, from commit `d91ed3d563019e649dc854c56ccd62551359488a` (`atac-seq/nucleosome-positioning`). See [`LICENSE`](LICENSE) for the preserved MIT notice.
