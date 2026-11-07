# Background

Not needed to run any route.

## An ASV is a denoiser's output on one run

The feature table is the residue of modeling decisions made before any result exists: which primers were stripped, where reads were truncated, what error model the run's quality scores supported, what was called a chimera. Declare these knobs. Illumina error rates are run-specific, so one model fit to concatenated runs denoises wrongly. One genome carries multiple, often divergent 16S copies (in about 60% of prokaryotes; *E. coli* K-12 has 7 copies in about 5 sequence types), so one organism can appear as several ASVs, while 97% OTUs lump the copies back.

## Tools

| Tool | Mechanism |
|------|-----------|
| DADA2 | Per-run parametric error model, abundance-partition denoising, merge, chimera; pseudo-pooling for singleton sensitivity |
| q2-dada2 | QIIME2 wrapper: denoise-paired, denoise-single, denoise-pyro, denoise-ccs |
| Deblur | Static upper-bound Illumina error profile (positive filter), one fixed length |
| cutadapt | Primer and adapter trimming |
| ITSxpress | HMM-trims the variable-length ITS spacer, keeping quality scores |
| VSEARCH | 97% OTU clustering, dereplication, chimera detection |

## Thresholds and their origin

| Threshold | Rationale |
|-----------|-----------|
| `maxEE` 2,2 (reverse to 5 for long amplicons) | Expected-errors filter beats a hard Q cutoff |
| truncLen sum >= amplicon length + 12 | `mergePairs` `minOverlap` default |
| truncLen <= post-primer read length | `filterAndTrim` drops, never pads, short reads |
| `maxN` 0 | DADA2 cannot model ambiguous bases |
| Chimera retained-read fraction 0.8-0.99 | Chimeras are many ASVs, few reads |
| `pool='pseudo'` | Approximates quadratic full pooling in linear time |
| Deblur one trim length | Positive filter needs a single read length |
| Copy-number correction | Predictable only near reference genomes; correction can add error |
