---
name: bio-cfdna-preprocessing
description: Decides how to preprocess plasma cfDNA sequencing data so the recoverable signal survives - library-prep-aware fragment expectations (dsDNA vs ssDNA/adaptase prep), UMI/duplex consensus with fgbio (ExtractUmisFromBam, GroupReadsByUmi --strategy paired for duplex, CallMolecularConsensusReads vs CallDuplexConsensusReads, FilterConsensusReads min-reads "total s1 s2"), the align->group->consensus->RE-align ordering, and the cfDNA dedup trap where naive coordinate dedup collapses nucleosome-coincident independent molecules. Covers assay-conditional simplex versus duplex choices, the singleton/sensitivity tax at low input, and reading the insert-size histogram as a pre-analytical QC instrument. Use when processing plasma cfDNA reads before fragmentomics, ctDNA mutation calling, or tumor-fraction estimation.
tool_type: mixed
primary_tool: fgbio
license: MIT
author: GPTomics
---

## Version compatibility

Reference examples target bwa 0.7.17+, fgbio 2.1+, numpy 1.26+, pysam
0.22+, and samtools 1.19+.

Before applying a recipe, inspect the installed versions and live help:

- Python: `pip show <package>` and `help(module.function)`
- CLI: `<tool> --version` and `<tool> --help`

fgbio flag semantics can drift. In particular, verify that
`CallDuplexConsensusReads --min-reads` is a permissive pre-filter and that
`FilterConsensusReads --min-reads` accepts `total strand1 strand2`, with the
more stringent strand value first when strand thresholds differ.

# cfDNA preprocessing

Preserve the physical cfDNA signal while producing analysis-ready reads.
Treat library preparation and pre-analytics as hard limits: consensus can
suppress errors, but it cannot recover molecules discarded by collection,
plasma preparation, or library chemistry.

## Establish the assay before choosing a path

Record:

1. library chemistry: dsDNA ligation, native-end ssDNA, or adaptase/tail-based
   ssDNA;
2. UMI design: none, single-strand, or duplex, including read structures and
   any skip/stem bases;
3. downstream aim: low-VAF mutation detection, MRD-grade specificity,
   fragmentomics, or sWGS/ULP-WGS tumor-fraction estimation;
4. input mass, observed family depths, and whether sensitivity or genotyping
   specificity is the priority;
5. whether physical or in-silico size selection has already conditioned the
   fragment distribution.

Use the detailed [method-selection and fragment-signal
reference](references/method-selection-and-fragment-qc.md) when chemistry,
consensus level, low-input recovery, or insert-size interpretation is in
question.

## Choose the processing path

| Scenario | Path | Non-negotiable reason |
|---|---|---|
| Targeted panel with single-strand UMIs | adjacency grouping -> molecular consensus -> consensus filtering | Suppresses PCR/sequencer error while retaining simplex families |
| MRD or other high-specificity low-VAF assay | Consider duplex prep -> paired grouping -> duplex consensus -> a validated two-strand filter | Duplex can reject template-resident damage, but use it only when empirical recovery, molecular depth, background error, target design, and the final calling rule support the assay's claimed LoD |
| sWGS/ULP-WGS tumor fraction | trim -> align -> light dedup; no consensus | Copy-number inference needs even coverage, not a consensus error floor |
| Fragmentomics or end motifs | native-end ssDNA prep; no upstream size selection | End repair and size selection alter the signal being measured |
| Picogram input, detection priority | call permissively and retain singleton evidence | Requiring duplicate observation can discard the only mutant molecule |
| No UMIs, quantitative cfDNA | do not coordinate-deduplicate by position alone | Nucleosome-positioned independent molecules often share endpoints |

## Run UMI consensus in the required order

For UMI-bearing targeted data, use this order:

1. Extract each inline UMI segment to a distinct molecular-index tag and retain
   the combined `RX` tag. Include the `S` skip/stem token where applicable
   (for example `6M11S+T`).
2. Query-group the uBAM, convert it to paired FASTQ, align with `bwa mem -Y -p`,
   and use `ZipperBams` to restore UMI and read metadata. Template-coordinate
   sort the zipped alignments for grouping.
3. Group by UMI and approximate position. Use `adjacency` for simplex and
   `paired` for duplex; paired grouping is mandatory for reconstructing the
   two duplex strands.
4. Call molecular or duplex consensus permissively. The caller emits unmapped
   consensus reads because the consensus sequence differs from every input
   read.
5. Re-align the consensus sequence and use `ZipperBams` to transfer consensus
   and UMI tags to the new alignments, then queryname-sort the mapped records.
6. Apply the real quality gate with `FilterConsensusReads` to query-grouped
   input, then coordinate-sort and index its output. A duplex
   `--min-reads 2 1 1` requires both strands; a permissive simplex path can use
   a single value or an explicitly documented strand policy.
7. Index the final BAM and preserve the family, consensus, chemistry, and
   filtering metadata needed downstream.

Use the complete command recipe and flag commentary in the [fgbio consensus
workflow](references/fgbio-consensus-workflow.md). Do not omit the second
alignment or substitute coordinate deduplication for UMI-family grouping.

## Inspect the fragment distribution before trusting downstream results

Summarize proper-pair, primary, PF, nonduplicate, positive template lengths. At
minimum report the observation count, mode, median, short fraction (90-150 bp),
fraction over 250 bp, and counts filtered for each flag or bound. Use
[`scripts/preprocess_cfdna.py`](scripts/preprocess_cfdna.py) for the reusable
`insert_size_qc()` implementation; it rejects a nonpositive `max_size`.

Interpret the distribution relative to library chemistry:

- a mode around 167 bp and a roughly 10.4 bp sub-mode periodicity are expected
  for a healthy mononucleosome cfDNA distribution;
- excess long-fragment mass and a mode shift upward can indicate leukocyte-gDNA
  contamination;
- an adaptase library can shift the mode roughly 10 bp lower without implying
  contamination;
- a 120-130 bp spike can be adapter dimer and should be addressed before
  interpreting the biological histogram.

Do not apply blanket high-MAPQ filtering or length selection before
fragmentomics: both preferentially delete the short, ctDNA-enriched tail.

## Runnable resources

- [`scripts/preprocess_cfdna.py`](scripts/preprocess_cfdna.py) exposes
  `preprocess_cfdna()` for the fgbio/bwa/samtools chain and
  `insert_size_qc()` for fragment-length QC. Direct execution prints an API
  summary; it is not a command-line argument parser.
- [fgbio consensus workflow](references/fgbio-consensus-workflow.md) contains
  the shell-level recipe for users assembling the pipeline manually.

The Python pipeline invokes `fgbio`, `bwa`, and `samtools` as connected argv
processes (never through a command shell) and imports `pysam` and `numpy`. It
expects a paired, unmapped UMI-bearing BAM, one molecular-index tag per `M`
segment, a compatible indexed reference FASTA, and writable space for
intermediate BAMs. Installation and prompt examples remain in the [usage
guide](usage-guide.md).

## Report and preserve

For every run, report the library chemistry, UMI/read structure, simplex versus
duplex choice, grouping strategy, caller and filter thresholds, input mass,
raw and consensus molecule counts, insert-size summary, reference build, tool
versions, and every size/MAPQ selection. Do not present fragment-derived
features from a length-selected library as unbiased.

## Routed detail

- [Method selection and fragment QC](references/method-selection-and-fragment-qc.md):
  chemistry landscape, consensus error floors, low-input tradeoffs, failure
  modes, and quantitative thresholds.
- [fgbio consensus workflow](references/fgbio-consensus-workflow.md): exact
  command sequence, flag semantics, and execution-specific errors.
- [Scientific references](references/scientific-references.md): source
  bibliography for the numerical and methodological claims.

## Related skills

- analytical-validation - LoD and molecule-counting framework fed by input QC
- fragment-analysis - fragmentomics consumes the preserved fragment ends
- ctdna-mutation-detection - consensus reads feed low-VAF calling
- tumor-fraction-estimation - sWGS minimal-processing path
- alignment-files/duplicate-handling - general deduplication context
- read-qc/quality-reports - upstream read QC
