# Method selection and fragment QC

Use this reference when choosing library chemistry, simplex versus duplex
consensus, family thresholds, or interpreting a cfDNA insert-size histogram.

## Pre-analytics and library chemistry set the ceiling

cfDNA is a nucleosome footprint, not randomly sheared DNA. Its fragment-length
distribution carries biological signal: a mononucleosome mode near 167 bp, a
short ctDNA-enriched tail around 134-144 bp, and sub-nucleosomal fragments.

Two upstream choices determine what survives:

1. Delayed blood processing can release leukocyte gDNA, dilute tumor fraction,
   and add long fragments. No downstream algorithm can recover the lost
   denominator. Route collection and plasma-preparation validation to the
   analytical-validation skill.
2. Library chemistry changes the recoverable distribution. dsDNA ligation
   requires duplex substrate and end repair, whereas native-end ssDNA methods
   retain short, nicked, damaged molecules. Adaptase chemistry can shift the
   apparent mode.

UMI consensus is error suppression, not just duplicate removal. A
single-strand consensus removes errors introduced after amplification begins,
but a lesion already present on the template (for example C-to-T deamination
or G-to-T oxidation) is copied into the entire family. Duplex consensus can
reject such damage because both original strands must agree.

## Chemistry and consensus landscape

| Choice | What it does | Recoverable distribution or error floor |
|---|---|---|
| dsDNA ligation prep (NEB/KAPA-style) | Requires duplex substrate; end repair and A-tailing polish ends | Clean roughly 167 bp mode; loses sub-100 bp and native-end signal |
| Native-end ssDNA prep (SRSLY/Kircher-Meyer lineage) | Denatures and ligates single strands | Recovers short, nicked, damaged, and sub-nucleosomal molecules |
| Adaptase/tail-based ssDNA prep (Swift/Accel-1S) | Uses single-strand adaptase chemistry | Mode can appear about 10 bp short; sawtooth may be blunted |
| Single-strand UMI consensus | Majority-votes reads from one source strand | Suppresses many PCR/sequencer errors; achieved background and LoD remain assay-specific |
| Duplex consensus | Combines strand-specific consensuses | Can reject strand-specific template damage; achieved background and LoD still depend on recovery, molecule count, targets, and calling rules |

## Scenario guide

| Scenario | Recommended path | Why |
|---|---|---|
| Deep targeted panel with single-strand UMIs | adjacency group -> molecular consensus -> consensus filter | Standard simplex error suppression for panel VAF work |
| High-specificity low-VAF or MRD assay | Evaluate duplex prep -> paired group -> duplex consensus -> empirically selected strand filter | Prefer duplex only when measured recovery, molecular depth, background suppression, target design, and the final caller improve the validated assay tradeoff |
| sWGS/ULP-WGS tumor fraction | trim -> align -> light dedup; no consensus | Copy-number signal depends on even coverage |
| Degraded or low-input sample | Compare validated dsDNA and ssDNA recovery for the specimen and target use | ssDNA can recover damaged or ultrashort molecules, but platform-specific background and conversion efficiency may change the best choice |
| Fragmentomics/end motifs | Prefer minimally perturbing chemistry and no upstream size selection after assay-specific validation | End repair, ligation chemistry, and selection can bias the features; no single chemistry is universally unbiased |
| Picogram input, detection priority | permissive family threshold; accept singletons | Duplicate-observation requirements impose a sensitivity tax |
| No UMIs, quantitative readout | avoid coordinate deduplication by position alone | Independent cfDNA molecules share nucleosome-positioned endpoints |

Methodology evolves. Confirm current fgbio semantics and vendor-specific
fragment distributions before freezing a pipeline. Never infer an achieved
assay LoD from a published molecular error floor alone: validate recovery,
background error, molecule counts, target design, pre-analytics, and the final
calling rule on representative samples.

## Insert-size QC

Use the reusable `insert_size_qc()` function in
[`../scripts/preprocess_cfdna.py`](../scripts/preprocess_cfdna.py). It counts
proper-pair, primary, PF, nonduplicate alignments with positive template
lengths up to a positive chosen maximum and reports:

- number of observations;
- modal and median fragment length;
- fraction in the 90-150 bp ctDNA-enriched window;
- fraction over 250 bp as a long-fragment indicator;
- overlapping counts for every excluded flag and length condition.

Interpret all cutoffs relative to chemistry and assay design:

- a mode near 167 bp plus about 10.4 bp sub-mode periodicity is the canonical
  mononucleosome pattern;
- excess long-fragment mass (roughly above 180-250 bp) and a mode shift upward
  can indicate leukocyte-gDNA contamination;
- a mode about 10 bp low can be expected for adaptase libraries rather than a
  correction target;
- a spike around 120-130 bp can indicate adapter dimer; trim adapters before
  reading the histogram biologically.

## Failure modes that alter scientific meaning

### Coordinate deduplication without UMIs

Independent cfDNA molecules cluster at nucleosome and linker boundaries, so
many legitimate molecules share start and end coordinates. Picard
`MarkDuplicates` or `samtools markdup` can label them as PCR duplicates,
deflating molecule counts and low-frequency VAF. Use UMI families when
available; without UMIs, avoid position-only deduplication for quantitative
readouts and disclose the bias.

### Strictness applied at the duplex caller

Raising `CallDuplexConsensusReads --min-reads` discards families at the
pre-filter, before `FilterConsensusReads` can apply the actual policy. Call
permissively and filter the re-aligned consensus.

### Adjacency grouping on duplex data

Adjacency grouping cannot pair the two complementary UMI orientations. The
apparent duplex output can collapse to simplex. Use paired grouping.

### Blanket high-MAPQ filtering

Short 40-80 bp fragments have fewer anchoring bases and are overrepresented in
the low-MAPQ tail even when correctly placed. A blanket MAPQ 30 or 60 rule can
preferentially delete the short, ctDNA-enriched molecules. Set and report a
threshold in the context of fragment length and downstream purpose.

### The singleton tax at low input

At picogram input, many unique molecules are observed only once. Requiring two
reads per family can discard the only evidence for a true low-VAF variant.
Favor molecular recovery for detection; reserve strict two-strand policies for
genotyping or MRD settings where specificity dominates.

### Fragmentomics after size selection

Physical or in-silico size selection conditions on the feature being measured.
It can enrich detection sensitivity but invalidates unbiased fragment-ratio,
end-motif, and length-distribution interpretation. Retain selection metadata
and do not report selected distributions as native.

## Quantitative anchors

| Anchor | Source | Interpretation |
|---|---|---|
| 167 bp mononucleosome mode; about 10.4 bp periodicity | Snyder 2016 | Nucleosome core plus linker and helical pitch |
| 35-80 bp transcription-factor/CTCF footprints | Snyder 2016 | Sub-nucleosomal protection, most visible with ssDNA prep |
| 134-144 bp ctDNA principal length versus 167 bp germline | Underhill 2016 | Tumor-derived fragments are shorter on average |
| 90-150 bp and 250-320 bp enrichment windows | Mouliere 2018 | Selection can enrich mutant allele fraction but biases fragmentomics |
| About 10.7x mitochondrial and 71.3x microbial cfDNA enrichment with ssDNA prep | Burnham 2016 | dsDNA prep misses much of the ultrashort fraction |
| About 3x barcode and 3x in-silico suppression, about 15x combined | Newman 2016 | Family consensus and background modeling are complementary |
| Duplex error floor below 1 per 1e7 nt | Kennedy 2014 | Two-strand concordance rejects strand-specific damage |
| Filter defaults: max read error 0.025, max base error 0.1, max no-calls 0.2 | fgbio documentation | Confirm fraction-versus-count and defaults in installed help |

See [scientific references](scientific-references.md) for full citations.
