---
name: bio-long-read-splicing
category: Data Analysis
description: "Analyzes alternative splicing from PacBio Iso-Seq (HiFi, Kinnex/MAS-Iso-seq) and Oxford Nanopore (direct cDNA, direct RNA, R10.4.1+) long-read RNA-seq with full-isoform resolution. Tools include FLAIR (correct/collapse/quantify/diffSplice for PacBio + ONT), IsoQuant (de-novo or annotation-guided isoform discovery 2024 SOTA), Bambu (annotation-aware Bayesian discovery + quantification with Novel Discovery Rate), SQANTI3 (isoform classification: FSM/ISM/NIC/NNC + artifact flags), rMATS-long (event calling on long-read isoforms), and minimap2 (-ax splice:hq for HiFi; -ax splice -k14 for ONT cDNA; add -uf only for direct RNA or stranded cDNA preps). Handles annotated microexons (with junction-guided alignment), recursive splicing, complex multi-exon isoforms, and DTU with lower (not zero) transcript-assignment uncertainty. Use when short-read AS limitations (anchor length, complex isoforms, microexons, recursive splicing, transcript ambiguity) demand full-isoform resolution."
tool_type: mixed
primary_tool: FLAIR
license: MIT
author: GPTomics
---

## Version Compatibility

Reference examples checked (2026-09) with: FLAIR 3.0.1, IsoQuant 4.0.0, Bambu 3.8.3 (loads its models only with xgboost 1.x), SQANTI3 6.0.2, minimap2 2.31, samtools 1.24, bedtools 2.31.1, gffread 0.12.9, rMATS-long 2.1.0, uLTRA 0.1, DRIMSeq 1.34.0, stageR 1.28.0. skera 1.4.0, lima 26.2.1, and isoseq 26.2.0 were run in sequence on a synthetic Kinnex array with valid `zm` tags (see the single-cell section); no real Kinnex sample was evaluated.

Flags changed between majors: FLAIR 2.x `flair correct` had `--genome` and `--shortread` (3.x: `--junction_tab`/`--junction_bed`, no genome); IsoQuant's entry point is `isoquant` (`isoquant.py` exists only in a git checkout); SQANTI3 6.x has no `--skipORF` (ORF prediction is off unless `--include_ORF`).

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Long-Read Splicing Analysis

Full-length long-read sequencing helps with problems that short-read AS cannot solve: anchor-length-limited microexon detection (with junction-guided alignment), complex multi-exon isoform deconvolution, recursive splicing in long introns, and transcript-quantification uncertainty in DTU. The 2024-2026 transition: long-read is becoming the splicing default for high-resolution analysis.

## Install

```bash
# bioconda: flair isoquant minimap2 samtools bedtools gffread sqanti3 rmats-long ultra_bioinformatics
conda install -c conda-forge -c bioconda flair isoquant minimap2 samtools bedtools gffread   # SQANTI3 (Python 3.11 + R) and rMATS-long: give each its own env
# R: BiocManager::install(c('bambu', 'IsoformSwitchAnalyzeR', 'DRIMSeq', 'stageR'))
# PacBio (single-cell block only): conda install -c bioconda pbskera lima isoseq
```

## When Long-Read Wins

| Question | Why long-read wins |
|----------|---------------------|
| Microexon detection (3-27 nt) | Reads span the microexon entirely, but minimap2 drops it unless the junctions are given (see Microexons) |
| Long-intron recursive splicing | Can detect ratchet point usage (Sibley 2015 *Nature*) |
| Complex isoform deconvolution (TTN, MAPT, NEFM) | Single read per isoform avoids EM ambiguity |
| DTU with lower quantification uncertainty | Reads are assigned to isoforms directly (no EM); truncated and ambiguous reads still leave some uncertainty (see DTU section) |
| Novel transcript discovery | No annotation dependence |
| Phasing splicing with SNVs | Allele-resolved isoforms |
| Single-cell full-length isoforms | MAS-Iso-seq + 10X 5' is the practical SOTA |
| Cryptic splicing in TDP-43 ALS | Full-length reads confirm cryptic exon inclusion in target transcripts |

## Analysis gates

- Before reporting SQANTI3 results, filter or explicitly review intra-priming and RT-switching flags; otherwise cDNA artifacts can look like novel isoforms or fusions.
- Start Bambu at `NDR = 0.1`; do not loosen it to chase novel calls without orthogonal validation. On sparse or few-gene inputs, treat its `NDR approximated` warning as a stop condition for novel-transcript claims.
- For a runnable end-to-end starting point, use [`examples/longread_splicing_pipeline.sh`](examples/longread_splicing_pipeline.sh); it takes platform and input paths through environment variables.

## Platform Selection Matrix

| Platform | Throughput | Accuracy (modal) | Best for | Fails when |
|----------|------------|------------------|----------|------------|
| PacBio Revio HiFi (Iso-Seq) | ~25M reads / SMRT cell | Q30+ (CCS) | Bulk transcript discovery; gold standard | Cost prohibitive for very large cohorts |
| PacBio Kinnex / MAS-Iso-seq | ~16x Iso-Seq via concatemer | Q30+ | High-throughput single-cell long-read | Kinnex de-array (skera) is an extra step |
| ONT direct cDNA (R10.4.1, PCS-114) | Millions / flowcell | ~98% simplex, ~99% duplex | Cost-effective; throughput | Minor higher error than HiFi |
| ONT direct RNA (RNA004, 2024+) | ~30M reads | ~96-98% | Native modifications (m6A, pseudo-U); no RT bias | Lower throughput; higher input |
| ONT pre-R10 (R9.4.1) | Same as R10 | ~85-90% | Legacy data | Pre-R10 not recommended for splicing analysis (false novel junctions) |

**Read length:** PacBio HiFi cdna typically 1-10 kb; ONT cdna 0.5-50+ kb (long-tailed). Both span typical mammalian transcripts. Direct RNA on ONT preserves true 5'/3' termini and modifications.

## Decision Tree by Use Case

| Use case | Recommended tools |
|----------|--------------------|
| Bulk Iso-Seq transcript discovery in well-annotated organism | minimap2 -ax splice:hq -> IsoQuant or Bambu -> SQANTI3 |
| Bulk ONT cDNA in well-annotated organism | minimap2 -ax splice -k14 (no -uf) -> IsoQuant or FLAIR -> SQANTI3 |
| Microexons (3-27 nt) | minimap2 --junc-bed (annotated; ONT also --junc-bonus 16) -> IsoQuant; junctions from short reads for unannotated ones (see Microexons) |
| End-to-end pipeline for differential analysis | FLAIR (correct -> collapse -> quantify -> diffSplice) |
| Joint discovery + quantification with calibrated novel rate | Bambu in R |
| De novo discovery for non-model organism | IsoQuant with --genedb omitted |
| Event-level differential splicing on long reads | rMATS-long |
| DTU on long-read transcript counts | FLAIR quantify -> DRIMSeq -> stageR (no Salmon Gibbs needed) |
| Hybrid short+long for cohort | FLAIR correct with short-read junctions (`--junction_tab`/`--junction_bed`) |
| Single-cell full-length isoforms | MAS-Iso-seq + 10X 5' -> skera -> lima -> isoseq refine (see single-cell-splicing) |
| Cryptic exon validation in ALS | minimap2 -> FLAIR collapse -> manual inspection of UNC13A, STMN2 |
| ASO design with full-isoform context | minimap2 -> IsoQuant -> SQANTI3 -> ASO design (see splice-variant-prediction) |

## Splice-Aware Alignment

```bash
# Annotation junctions for minimap2 (BED12); --junc-bed makes it prefer annotated junctions and rescues annotated microexons (ONT needs --junc-bonus 16 for the short ones, see Microexons)
gffread gencode.v45.annotation.gtf --bed -o annotation.bed12

# PacBio HiFi (Iso-Seq) -> minimap2 splice:hq preset
minimap2 -ax splice:hq --secondary=no --junc-bed annotation.bed12 \
    -t 16 \
    reference.fa \
    isoseq.fastq.gz | \
    samtools sort -@ 8 -o isoseq_aligned.bam
samtools index isoseq_aligned.bam

# ONT direct cDNA (PCS-114, PCB-114): reads come in both orientations, so no -uf
minimap2 -ax splice -k14 --secondary=no --junc-bed annotation.bed12 --junc-bonus 16 \
    -t 16 \
    reference.fa \
    ont_cdna.fastq.gz | \
    samtools sort -@ 8 -o ont_cdna_aligned.bam
samtools index ont_cdna_aligned.bam

# ONT direct RNA (RNA004): every read is in transcript orientation, so -uf is correct
minimap2 -ax splice -uf -k14 --secondary=no --junc-bed annotation.bed12 --junc-bonus 16 \
    -t 16 \
    reference.fa \
    ont_rna.fastq.gz | \
    samtools sort -@ 8 -o ont_rna_aligned.bam
samtools index ont_rna_aligned.bam
```

**`-uf` only for reads that are all in transcript orientation** (direct RNA, orientation-fixed FLNC, stranded preps). On mixed-orientation reads it forces half of them onto the wrong strand and creates false junctions. Measured on simulated reads (pysam intron chains vs planted truth): unstranded ONT cDNA `splice -uf -k14` 47% exact chains and 50% of reads with a false junction, versus 94% and 1.1% for `splice -k14`; unoriented HiFi with `splice:hq -uf`: 29% of reads with a false junction. On oriented reads (HiFi FLNC, direct RNA) `-uf` and no `-uf` gave identical chains, so omitting it costs nothing. On real LRGASP cDNA reads, `-uf` cut the reads on annotated junctions from 92% to 54-59% and FLAIR-corrected reads from 1536 to 972-1353.

Orientation check, on an alignment made **without** `-uf` (minimap2 tags spliced reads `ts:A:+` when the read is in transcript orientation):

```bash
samtools view -F 2308 aligned.bam | awk '{for(i=12;i<=NF;i++) if($i ~ /^ts:A:/){n++; if($i=="ts:A:+") p++}} END{printf "%.3f\n", p/n}'
# ~1.000 -> oriented reads (-uf is safe); ~0.5 -> unoriented (never -uf). Measured: oriented HiFi 1.000, unstranded ONT 0.501, real LRGASP cDNA 0.503
```

`--secondary=no` discards secondary alignments. `--junc-bed` (with or without `--junc-bonus 16`) does not block novel junctions: a novel 24-nt acceptor shift stayed correctly aligned. Without an annotation, drop `--junc-bed`.

`splice:hq` is the preset for HiFi; plain `splice -k14` is for ONT. The presets matter little on simulated reads (HiFi: same chains within 2 reads; ONT with `splice:hq`: 2.1% versus 1.1% of reads with a false junction), so the main alignment risk is `-uf`, not the preset.

### Microexons (3-27 nt)

Use annotation junctions with minimap2: **HiFi `--junc-bed` alone; ONT cDNA/direct RNA `--junc-bed --junc-bonus 16`**. A high bonus can recode reads that skip the exon into inclusion reads. Test both isoforms on a known locus before raising it: inclusion may rise, but skipping must not fall; intron-next-to-soft-clip reads must remain zero. This recipe was checked across 60 simulated 4-13 nt exons (200 inclusion + 200 skipping reads each) and independently on 12 sizes from 3-27 nt with tandem pairs. Direct-RNA results vary with error profile: at bonus 16, mean inclusion/skipping was 99.1/99.2% at ~4% simulated error and 93.2/94.9% at 7.2%; do not treat one error profile as a platform guarantee. On real LRGASP ONT cDNA, bonus 16 also changed 109/1883 (5.8%) chains by extending them through annotated introns, while keeping zero dangling junctions. See [microexon validation](references/microexons.md) for the measurements, controls, and short-read-junction route.

```bash
samtools view -F 2308 aligned.bam | awk '$6 ~ /N[0-9]+S$/ || $6 ~ /^[0-9]+S[0-9]+N/' | wc -l   # reads with an intron next to a soft clip
```

```python
# usage: python check.py aligned.bam UP_START UP_END DOWN_START DOWN_END SKIP_START SKIP_END
# intron coordinates 0-based half-open (BED12 block gaps): UP/DOWN = introns flanking the microexon, SKIP = exon1 end to exon3 start
import sys, pysam
bam, *c = sys.argv[1:]; c = list(map(int, c)); up, down, skip = tuple(c[0:2]), tuple(c[2:4]), tuple(c[4:6])
inc = exc = 0
for r in pysam.AlignmentFile(bam):
    if r.is_unmapped or r.is_secondary or r.is_supplementary: continue
    pos, introns = r.reference_start, set()
    for op, n in r.cigartuples:
        if op == 3: introns.add((pos, pos + n))
        if op in (0, 2, 3, 7, 8): pos += n
    inc += up in introns and down in introns
    exc += skip in introns
print("inclusion reads", inc, "skipping reads", exc)
```

Example output for a simulated 4-nt exon (200 inclusion + 200 skipping HiFi reads): plain `0 / 399`, default bonus `199 / 201`, bonus 20 `400 / 0` (every skipping read recoded).

An unannotated microexon needs externally supported junctions, for example a 6-column BED from STAR `SJ.out.tab`. Do not substitute uLTRA as an equivalent rescue: its 150/150 result was one annotated 10-nt sequence; on the independent 3-27 nt set it lost 3-nt exons and a 4+12 nt tandem pair and was weaker below about 12 nt, especially on ONT/direct RNA. See [microexon validation](references/microexons.md).

```bash
awk 'BEGIN{OFS="\t"} $4>0 {print $1,$2-1,$3,"sj"NR,$7,($4==1?"+":"-")}' SJ.out.tab > sr_junctions.bed
# HiFi; for ONT add --junc-bonus 16
minimap2 -ax splice:hq --secondary=no --junc-bed sr_junctions.bed -t 16 reference.fa isoseq.fastq.gz | samtools sort -o isoseq_aligned.bam
```

## FLAIR Workflow (correct -> collapse -> quantify -> diffSplice)

**Goal:** Identify, quantify, and test full-length isoforms from long-read RNA-seq across conditions.

**Approach:** Correct splice junctions against annotation and/or short-read junctions, collapse isoforms, quantify per-sample expression, run diffSplice for differential isoform usage.

```bash
# BAM -> BED12 (or `flair align`), one file per sample
bedtools bamtobed -bed12 -i aligned.bam > aligned.bed

# correct takes no genome; needs --gtf and/or short-read junctions (STAR SJ.out.tab via --junction_tab, or BED via --junction_bed)
flair correct \
    --query aligned.bed \
    --gtf gencode.v45.annotation.gtf \
    --junction_tab SJ.out.tab \
    --output flair_corrected \
    --threads 16

# collapse all samples together: cat the corrected BEDs, pass every reads file
flair collapse \
    --query flair_corrected_all_corrected.bed \
    --reads sample.fastq.gz \
    --genome reference.fa \
    --gtf gencode.v45.annotation.gtf \
    --output flair_collapsed \
    --threads 16

# reads_manifest.tsv: no header, tab-separated: sample_id  condition  batch  /path/reads.fastq.gz
flair quantify \
    --reads_manifest reads_manifest.tsv \
    --isoforms flair_collapsed.isoforms.fa \
    --output flair_quantified \
    --threads 16

flair diffSplice \
    --isoforms flair_collapsed.isoforms.bed \
    --counts_matrix flair_quantified.counts.tsv \
    --out_dir flair_diffsplice \
    --test \
    --threads 16
```

FLAIR (Tang 2020 *Nat Commun*) handles ONT and PacBio with the same workflow. `flair quantify` writes `<output>.counts.tsv` (first column `ids` = `<isoform>_<gene>`, sample columns `<id>_<condition>_<batch>`; `--sample_id_only` gives `ID` and plain sample ids). `flair diffSplice` writes per-event-type (`es`, `alt5`, `alt3`, `ir`) inclusion/exclusion tables; `--out_dir` must not exist yet (`-of` overwrites). `--test` adds a DRIMSeq test per event type and needs an `Rscript` on PATH with DRIMSeq, argparse, and data.table; pip/conda FLAIR provides none of them. An event type may have rows but fail DRIMSeq filtering with `No genes left after filtering!`; FLAIR can still exit 0. For every event type required by the analysis, require its `drimseq_<event>_*.tsv` output rather than treating process exit status as success.

## IsoQuant for Discovery + Quantification

**Goal:** De novo or annotation-guided isoform discovery and quantification with high precision.

**Approach:** Run `isoquant` (not `isoquant.py`) with reference + reads (or a minimap2 BAM via `--bam`) + data type; output is GTF + counts.

```bash
isoquant \
    --reference reference.fa \
    --genedb gencode.v45.annotation.gtf \
    --fastq sample1.fastq.gz sample2.fastq.gz \
    --data_type pacbio_ccs \
    --output isoquant_output \
    --threads 16 \
    --model_construction_strategy default_pacbio
```

`--data_type` accepts `pacbio_ccs` (HiFi), `pacbio`, `nanopore`/`ont`, `assembly` or `transcripts`; `--stranded forward` for reads in transcript orientation (direct RNA), default `none`. `--genedb` is optional for de novo discovery. Outputs are in `<output>/<prefix>/`: `<prefix>.transcript_models.gtf` (also holds gene, CDS, UTR and codon rows copied from the annotation), `<prefix>.transcript_counts.tsv` (annotated transcripts only, pooled over all input files; default prefix `OUT`) and `<prefix>.discovered_transcript_counts.tsv` (annotated transcripts plus the novel models, ids like `transcript115.chrQ.nnic`), each with a `..._grouped_file_name_counts.tsv` twin (one column per input file; its id column is headed `gene_id` although the rows are transcripts). 4.0.0 writes no `transcript_model_counts.tsv`. Transcript counts default to `--transcript_quantification unique_only`, so ambiguous reads are not counted (`__ambiguous` row: 22 of 300 reads in the ONT microexon test). A novel splice site shifted by only ~24 nt from an annotated one is merged into the annotated isoform (all 20 planted reads of one sample counted in it, +13%); confirm such sites with short-read junctions (FLAIR `--junction_tab`). A planted 30-nt donor shift had no model in a `--genedb` run either, and `--illumina_bam` (125,055 simulated short reads, 40,278 spliced) changed nothing; a run without `--genedb` modelled it (43 reads), as it did an unannotated exon skip (37 reads, also modelled with `--genedb`). IsoQuant (Prjibelski 2023 *Nat Biotech*) is current SOTA for novel transcript reconstruction; pairs well with SQANTI3 for downstream classification.

Memory requirement: >=64 GB for atlas-scale runs.

## Bambu for Annotation-Aware Discovery + Quantification

**Goal:** Joint discovery and quantification with statistical filtering of novel isoforms.

**Approach:** R Bioconductor package; takes BAM + reference annotation + genome; outputs ranged SE objects of known + novel transcripts.

```r
library(bambu)

bam_files <- c('sample1.bam', 'sample2.bam', 'sample3.bam')
genome <- 'reference.fa'
gtf <- 'gencode.v45.annotation.gtf'

bambuAnnotations <- prepareAnnotations(gtf)

se <- bambu(
    reads = bam_files,
    annotations = bambuAnnotations,
    genome = genome,
    NDR = 0.1,
    ncore = 8   # BiocParallel workers; on Windows R use ncore = 1
)

writeBambuOutput(se, path = 'bambu_output/')

tx_counts <- as.data.frame(assays(se)$counts)
gene_counts <- transcriptToGeneExpression(se)
```

Bambu (Chen 2023 *Nat Methods* 20:1187-1195) uses **NDR** (Novel Discovery Rate) as a single, calibrated parameter replacing per-sample heuristics:

| NDR | Interpretation |
|-----|----------------|
| 0.05 | Stringent; few novel transcripts; highest precision |
| 0.1 | Balanced (default) |
| 0.2-0.3 | Permissive; more novel discoveries; recall over precision |

Excellent for combined discovery + quantification when statistical filtering matters. With several BAMs on Windows R, `ncore = 8` stopped with "BiocParallel errors ... could not find function seqlengths" and `ncore = 1` ran (single-BAM runs worked with either; multicore was not tried on Linux). NDR needs enough read classes: Bambu warns "NDR approximated" below ~50 read classes, and on a 3.4k-read simulated set NDR 0.05 and 0.1 reported no novel transcript (0.3: one; 1.0: two). Use it on genome-scale data, not a few genes.

## SQANTI3 Classification

**Goal:** Classify discovered isoforms relative to reference; flag artifacts (intra-priming, RT-switching).

**Approach:** Run `sqanti3_qc.py` on the isoform GTF; review classification (FSM/ISM/NIC/NNC/antisense/genic/intergenic/fusion) and quality flags.

```bash
# isoforms.gtf: transcript and exon rows only (FLAIR isoforms.gtf as is; from IsoQuant:
#   awk -F'\t' '$3=="transcript" || $3=="exon"' <prefix>.transcript_models.gtf > isoforms.gtf)
# -o is a file prefix, -d the output directory. --CAGE_peak and --polyA_motif_list are optional: drop them without files
sqanti3_qc.py \
    --isoforms isoforms.gtf \
    --refGTF gencode.v45.annotation.gtf \
    --refFasta reference.fa \
    -o sqanti3 -d sqanti3_qc \
    --aligner_choice minimap2 \
    --CAGE_peak hg38.cage_peak_phase1and2combined_coord.bed \
    --polyA_motif_list human.polyA.list.txt \
    --cpus 8
# classification: sqanti3_qc/sqanti3_classification.txt

sqanti3_filter.py rules \
    --sqanti_class sqanti3_qc/sqanti3_classification.txt \
    --filter_gtf isoforms.gtf \
    -o sqanti3_filtered -d sqanti3_filtered \
    --skip_report
# filtered GTF: sqanti3_filtered/sqanti3_filtered.filtered.gtf; per-isoform filter_result: *_RulesFilter_classification.txt
```

Optional inputs: CAGE peaks (hg38) from refTSS (reftss.riken.jp) or the FANTOM5 file the SQANTI3 docs link (`Magdoll/images_public`, `SQANTI2_support_data/hg38.cage_peak_phase1and2combined_coord.bed.gz`, gunzip first); the polyA motif list from the same repo (`SQANTI2_support_data/human.polyA.list.txt`) or `data/polyA_motifs/` of the SQANTI3 GitHub repo (the conda package ships neither). Both ran on real hg38 isoforms (`within_CAGE_peak` and `polyA_motif_found` filled). `--report skip` on the QC step avoids the HTML/PDF report (the default report also ran, 33 s). `--skip_report` on the filter step is needed on small inputs: its R report failed (tidyselect "subscript out of bounds", exit 1) on 6 and 10 isoforms after the filter tables were written. SQANTI3 6.x has no `--skipORF`.

| SQANTI category | Meaning |
|-----------------|---------|
| FSM (Full Splice Match) | All junctions match reference |
| ISM (Incomplete Splice Match) | Subset of reference junctions |
| NIC (Novel In Catalog) | Novel combination of known junctions |
| NNC (Novel Not in Catalog) | Contains novel junction |
| Antisense | Overlaps gene on opposite strand |
| Genic | Within gene but no junction match |
| Genic intron | Entirely within an intron of a reference gene |
| Intergenic | Between genes |
| Fusion | Spans multiple genes |

SQANTI3 (Pardo-Palacios 2024 *Nat Methods* 21:793-797) is the long-read isoform-curation/QC tool, with structural categories and QC tailored to ONT/PacBio error patterns. **Filter intra-priming and RT-switching** flags before reporting.

## rMATS-long for Differential Isoform Analysis on Long-Read Data

**Goal:** Apply differential isoform analysis to long-read transcript abundance with classification and visualization.

**Approach:** rMATS-long is a multi-script Python pipeline distributed via bioconda; entry point is `rmats-long` followed by the script name. It supports two modes: **abundance-based** (using ESPRESSO-style abundance estimates) and **ASM-based** (Alternative Splicing Modules — sets of isoforms sharing exon-junction structure). Run preprocessing scripts in order before `rmats_long.py`.

```bash
# conda install -c conda-forge -c bioconda rmats-long   (own env)
set -euo pipefail   # the steps below otherwise carry on after a failure and write header-only tables
mkdir -p alignment_info

# Preprocessing pipeline (ASM mode); per-script flag names verified vs Xinglab/rmats-long
rmats-long organize_gene_info_by_chr.py --gtf annotation.gtf --out-dir gene_info_by_chr/

# simplify_alignment_info processes one sorted, indexed BAM at a time -> one TSV;
# samples.tsv (sample_id<TAB>tsv_path) is what organize_alignment_info_by_gene_and_chr.py reads
: > samples.tsv
for bam in *.bam; do
    id="${bam%.bam}"
    rmats-long simplify_alignment_info.py --in-file "$bam" --out-tsv "alignment_info/${id}.tsv"
    printf '%s\talignment_info/%s.tsv\n' "$id" "$id" >> samples.tsv
done
# every per-sample table must exist and be non-empty
while IFS=$'\t' read -r id tsv; do [ -s "$tsv" ] || { echo "empty or missing $tsv" >&2; exit 1; }; done < samples.tsv

rmats-long organize_alignment_info_by_gene_and_chr.py \
    --gtf-dir gene_info_by_chr/ \
    --out-dir organized/ \
    --samples-tsv samples.tsv

rmats-long detect_splicing_events.py --align-dir organized/ --gtf-dir gene_info_by_chr/ --out-dir events/
rmats-long create_gtf_from_asm_definitions.py --event-dir events/ --out-gtf asm.gtf
rmats-long count_reads_for_asms.py --align-dir organized/ --event-dir events/ --gtf-dir gene_info_by_chr/ --out-dir asm_counts/

# Main differential analysis (ASM mode)
# --group-1 / --group-2 each take the PATH to a file whose single line is a
# comma-separated list of sample IDs (matching the BAM basenames in --align-dir).
echo 'ctrl1,ctrl2,ctrl3' > group1.txt
echo 'trt1,trt2,trt3' > group2.txt
rmats-long rmats_long.py \
    --group-1 group1.txt \
    --group-2 group2.txt \
    --event-dir events/ \
    --asm-counts-dir asm_counts/ \
    --align-dir organized/ \
    --gtf-dir gene_info_by_chr/ \
    --out-dir rmats_long_output/ \
    --adj-pvalue 0.05 \
    --delta-proportion 0.05 \
    --average-reads-per-group 10

# the result tables must have data rows, not just a header
for t in differential_asms.tsv differential_isoforms.tsv; do
    [ "$(wc -l < "rmats_long_output/$t")" -gt 1 ] || { echo "rmats_long_output/$t has no rows" >&2; exit 1; }
done
```

Alternative, abundance-based mode (when you already have ESPRESSO-style estimates), instead of the ASM steps:

```bash
rmats-long rmats_long.py --abundance abundance.esp --updated-gtf updated.gtf \
    --group-1 group1.txt --group-2 group2.txt --out-dir rmats_long_output/ --no-splice-graph-plot
```

Key flags: `--adj-pvalue` (default 0.05), `--delta-proportion` (default 0.05), `--average-reads-per-group` (default 10), `--no-splice-graph-plot` (skip expensive splice-graph rendering).

rMATS-long is a separate tool from short-read rMATS-turbo. The predecessor `lr2rmats` used long reads only to *augment* the short-read rMATS GTF. The ASM framework treats AS as a set-of-isoforms problem, more natural for long-read data than rMATS-turbo's pre-defined event categories.

## DTU on Long-Read Counts

**Goal:** Apply a DRIMSeq + stageR DTU pipeline to long-read transcript counts (for the DEXSeq or satuRn alternatives see `isoform-switching`).

**Approach:** Use FLAIR (or Bambu) transcript counts as input; reads are assigned to isoforms directly, so no Salmon Gibbs samples or EM are needed. Assignment uncertainty is lower than with short reads but not zero: 5'-truncated and ambiguous reads remain (IsoQuant's default `unique_only` leaves some reads out of the counts; a novel splice site can be absorbed by its annotated neighbour, see IsoQuant). Ran on real FLAIR `counts.tsv` output (6 samples) and on a 300-gene planted DTU set (30/30 planted genes recovered, 6 extra genes called).

```r
library(DRIMSeq); library(stageR)   # do not attach DEXSeq here: it masks results()

manifest <- read.table('reads_manifest.tsv', sep = '\t', col.names = c('sample_id', 'condition', 'batch', 'reads'))
counts <- read.table('flair_quantified.counts.tsv', header = TRUE, sep = '\t', check.names = FALSE)

# First column = <isoform>_<gene> ('ids', or 'ID' with --sample_id_only); sample columns are <id>_<condition>_<batch> (just <id> with --sample_id_only)
sample_cols <- vapply(seq_len(nrow(manifest)), function(i) {
    hit <- intersect(c(manifest$sample_id[i], paste(manifest$sample_id[i], manifest$condition[i], manifest$batch[i], sep = '_')),
                     colnames(counts))
    stopifnot(length(hit) == 1)
    hit
}, character(1))
dm_counts <- data.frame(gene_id = sub('^.*_', '', counts[[1]]), feature_id = counts[[1]],
                        setNames(counts[, sample_cols], manifest$sample_id), check.names = FALSE)
samples <- data.frame(sample_id = manifest$sample_id, condition = factor(manifest$condition))

n <- nrow(samples); n_min <- min(table(samples$condition))
d <- dmDSdata(counts = dm_counts, samples = samples)
d <- dmFilter(d, min_samps_feature_expr = n_min, min_feature_expr = 5,
              min_samps_feature_prop = n_min, min_feature_prop = 0.1,
              min_samps_gene_expr = n, min_gene_expr = 10)

design <- model.matrix(~ condition, data = DRIMSeq::samples(d))
d <- dmPrecision(d, design = design)
d <- dmFit(d, design = design)
d <- dmTest(d, coef = colnames(design)[2])
res_gene <- DRIMSeq::results(d)
res_tx <- DRIMSeq::results(d, level = 'feature')

# stageR: gene-level screen, transcript-level confirmation
res_gene <- res_gene[!is.na(res_gene$pvalue), ]
res_tx <- res_tx[res_tx$gene_id %in% res_gene$gene_id & !is.na(res_tx$pvalue), ]
pScreen <- setNames(res_gene$pvalue, res_gene$gene_id)
pConfirm <- matrix(res_tx$pvalue, ncol = 1, dimnames = list(res_tx$feature_id, 'transcript'))
tx2gene <- data.frame(transcript = res_tx$feature_id, gene = res_tx$gene_id)
sr <- stageRTx(pScreen = pScreen, pConfirmation = pConfirm, pScreenAdjusted = FALSE, tx2gene = tx2gene)
sr <- stageWiseAdjustment(sr, method = 'dtu', alpha = 0.05, allowNA = TRUE)
dtu <- getAdjustedPValues(sr, order = TRUE, onlySignificantGenes = TRUE)   # stage-wise adjusted p per gene and transcript; NULL when no gene is significant
```

The gene id is the text after the last `_` of the FLAIR id, which assumes gene ids without underscores. IsoformSwitchAnalyzeR v2 has explicit long-read input support.

## Single-Cell Long-Read for Splicing

**Goal:** Combine cell typing (10X 5' short read) with full-length isoform structure (Kinnex / MAS-Iso-seq).

**Approach:** Split 10X library; sequence half short-read for cell typing, half PacBio Kinnex for isoforms; de-array the Kinnex reads with skera, then run the Iso-Seq steps (lima, isoseq refine, isoseq cluster2). Barcode assignment to cells is outside this Skill (see `single-cell-splicing`).

```bash
# De-array Kinnex/MAS HiFi reads into segmented reads (S-reads). Obtain the
# kit-specific adapter FASTA from PacBio for the exact kit; do not guess its name.
skera split \
    raw_kinnex.bam \
    kit_adapters.fasta \
    segmented.bam

# Then lima -> isoseq refine -> isoseq cluster2 (bioconda lima, isoseq; see long-read-sequencing/isoseq-analysis)
```

skera 1.4.0 was run on a synthetic array. With a valid `zm` tag on its toy BAM, `skera -> lima --isoseq -> isoseq refine` completed (132 input S-reads -> 132 lima -> 132 FLNC). Without `zm`, the toy read names used `?` and lima stalled; this was a fixture defect, not a demonstrated lima limitation. No real Kinnex data was run.

Joglekar et al 2024 (*Nat Neurosci* 27:1051-1063) used this approach to map single-cell isoforms across developing and adult mouse and human brain. See `single-cell-splicing` for tools that work on the demultiplexed data.

## Per-Tool Failure Modes

### minimap2: -uf on Unoriented Reads

**Trigger:** `-uf` on ONT cDNA or unoriented HiFi (the recipes above give the rule).

**Symptom:** ~half the reads with false junctions, inferred gene strand wrong for ~half, few reads kept by `flair correct`; the orientation check gives ~0.5.

**Fix:** Drop `-uf` (the presets `splice:hq` and `splice -k14` differ little; see the alignment section).

### IsoQuant: Memory Pressure

**Trigger:** Atlas-scale cohort or low-RAM environment.

**Mechanism:** IsoQuant builds graph structures across all reads simultaneously.

**Symptom:** OOM kill; very slow runtime.

**Fix:** Increase RAM to >=64 GB; or batch by chromosome.

### Bambu: NDR Mistuning

**Trigger:** NDR=0.5+ or NDR=0.01.

**Mechanism:** NDR controls the precision-recall tradeoff for novel transcripts.

**Symptom:** Too many spurious novel transcripts (high NDR) or missing real novel transcripts (low NDR).

**Fix:** Default NDR=0.1 is balanced; adjust based on validation expectations.

### SQANTI3: RT-Switching Flags

**Trigger:** PacBio/ONT cDNA libraries with template switching artifacts.

**Mechanism:** RT-switching produces chimeric reads spanning two unrelated transcripts; SQANTI3 flags these.

**Symptom:** Many "fusion" transcripts in non-cancer samples; biologically implausible.

**Fix:** Filter out RT-switching flags via `sqanti3_filter.py`; investigate library prep if rate >5%.

### FLAIR: Novel Splice Sites Need Orthogonal Junctions

**Trigger:** Running `flair correct` with only `--gtf` (no `--junction_tab`/`--junction_bed`).

**Mechanism:** FLAIR keeps only junctions supported by the annotation or by the orthogonal (short-read) junctions; a novel splice site is moved to the nearest supported one or the read is dropped as inconsistent.

**Symptom:** A planted 24-nt novel-acceptor isoform (20 reads in one sample) was missing from `flair collapse`, and the annotated isoform's count was inflated by exactly those reads (194 vs 174).

**Fix:** Add the short-read junctions (STAR `SJ.out.tab` via `--junction_tab`, or a BED via `--junction_bed`, a 6-column BED: chrom, start, end, name, score = read support, strand). With them all 7 planted isoforms were recovered with exact counts.

## Reconciliation: When Long-Read Tools Disagree

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| FLAIR has more isoforms than IsoQuant | FLAIR collapse less stringent; or IsoQuant filtered more aggressively | Both tools have valid pipelines; report based on use case |
| Bambu calls fewer novel than IsoQuant | Bambu NDR=0.1 is more conservative | Adjust NDR or trust Bambu's calibration |
| SQANTI3 classifies as NNC, FLAIR thinks FSM | GENCODE version mismatch | Verify both tools use same annotation |
| Long-read isoform calls don't match short-read events | Short-read EM ambiguity; or long-read coverage gap | Trust long-read for unambiguous; trust short-read for high-coverage events |

## Quality Control for Long-Read Splicing

| Metric | PacBio HiFi | ONT cDNA R10.4.1 |
|--------|-------------|-------------------|
| Read accuracy (modal) | Q30+ (>=99.9%) | ~98% simplex / ~99% duplex |
| Splice junction concordance to short-read truth (unsourced estimate) | ~98% | 95-98% |
| Median read length (transcripts) | 1-4 kb | 0.5-3 kb |
| Throughput per run | ~25M HiFi reads | Tens of millions |
| Library input | 100-500 ng total RNA | 100-500 ng |
| Read direction | TSO + dT primed | TSO or random hexamer |

Figures are approximate vendor or literature values, not measured here.

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| `Bambu: Input annotation file not readable. Requires .gtf/.gff format or TxDb object` (`prepareAnnotations`) | Empty GTF or rows with fewer than 9 tab-separated columns (reproduced with rows cut to 6 columns; `gffread -E` did not flag them) | `awk -F'\t' '!/^#/ && NF!=9' annotation.gtf` lists the bad rows |

## Quality Thresholds

| Metric | Recommendation | Source |
|--------|----------------|--------|
| Full-length non-chimeric (FLNC) % | >=80% (PacBio Iso-Seq) | PacBio convention |
| FSM% | >=50% in well-annotated genome (field-convention rule of thumb; not specified in the SQANTI paper; real FLAIR-test isoforms gave 7 of 36 FSM, so it depends on annotation and sample) | SQANTI3 documentation; Tardaguila 2018 *Genome Res* 28:396 |
| NNC% | <=30% (>30% suggests artifacts unless biologically interesting) | SQANTI3 convention |
| Junction support | >=2 reads (or >=3 with strict filtering) | Conservative |
| Bambu NDR | 0.1 default; 0.05 stringent | Chen 2023 *Nat Methods* 20:1187 |
| SQANTI3 RT-switching flag | filter out unless validated | SQANTI3 convention |
| SQANTI3 intra-priming flag | filter out | SQANTI3 convention |
| HiFi CCS passes | >=3 | PacBio convention for Q30+ |

## Common Pitfalls

- **Ignoring reference annotation completeness** — SQANTI3 NNC categorization differs by GENCODE version; report version with results.
- **Not running `isoseq refine`** — concatemers and polyA artifacts inflate isoform counts.
- **Confusing FLAIR's 'collapse' with 'cluster'** — collapse merges similar isoforms post-alignment; `isoseq cluster2` merges raw reads pre-alignment.
- **Skipping CAGE / polyA validation in SQANTI3** — TSS / TTS hallucination is common in long-read isoforms.
- **DTU on too few replicates** — long-read is expensive; n=2 vs n=2 is common but underpowered.

## Related Skills

- splicing-quantification - Short-read PSI for cross-validation
- isoform-switching - DTU framework on long-read counts
- single-cell-splicing - MAS-Iso-seq + 10X integration
- long-read-sequencing/isoseq-analysis - PacBio Iso-Seq general pipeline (CCS, lima, refine, cluster)
- long-read-sequencing/long-read-alignment - minimap2 splice:hq details
- long-read-sequencing/long-read-qc - QC for long-read data
- splice-variant-prediction - Cross-reference variant predictions with full isoforms

## References

- Tang et al 2020 *Nat Commun* - FLAIR
- Prjibelski et al 2023 *Nat Biotech* - IsoQuant
- Chen et al 2023 *Nat Methods* 20:1187-1195 - Bambu
- Tardaguila et al 2018 *Genome Res* - SQANTI (original)
- Pardo-Palacios et al 2024 *Nat Methods* 21:793-797 - SQANTI3
- Pardo-Palacios et al 2024 *Nat Methods* 21:1349-1363 - LRGASP benchmark
- Wyman et al 2020 *bioRxiv* - TALON (note: not formally peer-reviewed)
- Li 2018 / 2021 *Bioinformatics* - minimap2
- Sahlin & Makinen 2021 *Bioinformatics* - uLTRA
- Sibley et al 2015 *Nature* - recursive splicing
- Al'Khafaji et al 2024 *Nat Biotech* - MAS-Iso-seq / Kinnex
- Joglekar et al 2024 *Nat Neurosci* 27:1051-1063 - scISOr-Seq2 single-cell brain isoform mapping
- Brown et al 2022 *Nature* - UNC13A cryptic exon (TDP-43)
- Klim et al 2019 *Nat Neurosci* - STMN2 cryptic splicing
