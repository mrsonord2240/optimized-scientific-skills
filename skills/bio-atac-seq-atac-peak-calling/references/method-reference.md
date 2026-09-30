# ATAC-seq Peak Calling - Method Reference

## Algorithmic Taxonomy

| Tool | Model | Treats fragments as | Min reps | Strength | Fails when |
|------|-------|---------------------|----------|----------|------------|
| MACS3/MACS2 | Local Poisson lambda + FDR | Point-source insertions (+/- shift) | 1 | Mature, ENCODE-default, fast, narrow + broad modes | Confounds NFR with broad accessible domains; no input means lambda from local genome only |
| Genrich (ATAC mode -j) | q-value on log-transformed p-value, joint replicate model (needs name-sorted BAMs) | Whole fragments (paired-end intervals) | 1 (multi-rep optional) | Treats reps jointly; can exclude chrM via `-e chrM`; auto blacklist via `-E`; PCR-dup removal via `-r` | Less peer-reviewed than MACS; thin literature; slow on deep libraries |
| MACS3 hmmratac (was HMMRATAC) | 3-state HMM (open / nucleosomal / background) on fragment-size signal | Fragment-size classes | 1 | Models nucleosome periodicity directly; differentiates NFR and flanking nucleosomes | Needs deep libraries (see failure modes); memory-hungry; slow; flat fragment distribution -> garbage HMM |
| HOMER `findPeaks -style dnase` | Fixed window + fold-change cutoff | Tag positions | 1 | Convenient for downstream HOMER motif analysis | Less calibrated p-values than MACS; window-size sensitive |
| nf-core/atacseq | Wrapper (MACS2 by default) | Same as MACS2 | 1 | Reproducible Nextflow pipeline with QC built in | Only as good as the underlying caller |

Methodology evolves; verify the current ENCODE ATAC-seq Standards (encodeproject.org pipelines/atac-seq) before locking parameters. ENCODE 4 still defaults to MACS2 (not MACS3) at time of writing; `macs3 callpeak` is API-compatible for ATAC parameters but not yet the official ENCODE binary.

## Shift-Extend vs BAMPE: The Critical Choice

Two valid ways to feed paired-end ATAC into MACS:

**Pattern A (ENCODE / "single-end-ified"):** `-f BAMPE` actually IGNORES `--shift/--extsize`. To activate them, use `-f BAM` and treat each end independently. ENCODE's pipeline uses `-f BAM --shift -75 --extsize 150` to model each Tn5 cut as a 150 bp window centered on the insertion site, ignoring fragment lengths.

**Pattern B (paired-fragment):** `-f BAMPE` uses the full paired-end fragment span as the signal interval. Best when fragment lengths are biologically meaningful (e.g., NFR-only peak calling at 38/75 bp). In BAMPE mode, do NOT set `--shift/--extsize` (silently ignored, but confusing).

For most bulk ATAC, Pattern A matches ENCODE convention and is reproducible against published peak sets. Pattern B can be more sensitive at narrow regulatory elements but does not match ENCODE outputs.

## Effective Genome Size

`-g hs` and `-g mm` are MACS shorthands for old defaults. Modern values:

| Genome | MACS shorthand | Actual mappable size | Source |
|--------|---------------|----------------------|--------|
| hg38 | `-g hs` (2.7e9) | 2.701e9 (50bp), 2.748e9 (75bp), 2.806e9 (100bp), 2.862e9 (150bp) | deepTools `effectiveGenomeSize` |
| hg19 | `-g hs` (2.7e9) | 2.686e9 (50bp), 2.777e9 (100bp) | deepTools |
| mm10 | `-g mm` (1.87e9) | 2.308e9 (50bp), 2.408e9 (75bp), 2.467e9 (100bp) | deepTools |
| mm39 | none | 2.310e9 (50bp), 2.468e9 (100bp) | deepTools |

Wrong size shifts every q-value but rarely changes peak ranks. Use `unique-kmers.py` (khmer) or the deepTools tabulated values for exact sizes; the shorthand is a decade-old approximation.

### When It Matters

**Trigger:** Comparing peaks across genome builds or species; reproducing published q-value cutoffs; hi-resolution lambda estimation.

**Mechanism:** MACS estimates genome-wide lambda as `total_reads / effective_size`. Wrong size -> wrong null -> shifted q-values, especially at the marginal cutoff.

**Symptom:** Peak counts diverge ~10-20% from published numbers when re-running an old dataset.

**Fix:** Pull the read-length-matched value from deepTools `effectiveGenomeSize` table. For pipelines, parameterize this; never inline the shorthand for cross-study comparisons.

## Per-Tool Failure Modes

### MACS2/MACS3 -- Confounded NFR + broad accessibility

**Trigger:** Cell type with extended open domains (e.g., active super-enhancers, MYOD1 regulons, locus-control regions).

**Mechanism:** Default narrow-peak mode segments wide accessible domains into multiple smaller peaks at local lambda spikes; `--broad --broad-cutoff 0.1` merges them but inflates total length and breaks IDR comparability.

**Symptom:** Peak count >> 200k for human bulk ATAC at ENCODE depth; mean peak width < 200 bp; visual inspection in IGV shows 3-5 calls under one continuous accessibility block.

**Fix:** Run both narrow and broad; use narrow for differential analysis, broad for domain-level enrichment (e.g., super-enhancer overlap). Do NOT use `--call-summits` for broad mode.

### Genrich -- Replicate weighting and chrM exclusion

**Trigger:** Replicates with very different library sizes; high-mitochondrial samples not pre-filtered.

**Mechanism:** Genrich's joint mode combines the per-replicate p-values at each position via Fisher's method. Library-size imbalance dominates the joint p-value; chrM reads inflate background unless `-e chrM` is set.

**Symptom:** Most-significant peaks cluster on chrM or on the largest-library replicate's high-coverage regions.

**Fix:** Always pass `-e chrM` (Genrich 0.6+) and `-E blacklist.bed`. Down-sample BAMs to common depth (`samtools view -s`) before joint calling if libraries differ >2x. Add `-r` to remove PCR duplicates inside Genrich, OR pre-deduplicate (do not do both). Inputs must be name-sorted (`samtools sort -n`); coordinate-sorted BAMs abort with `not sorted by queryname`. A tested command is in `usage-guide.md`.

### MACS3 hmmratac (HMMRATAC) -- Depth and fragment-size dependence

**Trigger:** Library < 25M nuclear reads, or libraries with degraded chromatin and flat fragment-size distribution.

**Mechanism:** The 3-state HMM is trained from fragment-size classes (NFR ~50 bp, mono ~200 bp, di ~400 bp peaks). Without periodicity the emission distributions collapse and the HMM cannot separate states.

**Symptom:** The `_accessible_regions.narrowPeak` output is empty, or all regions are tiny (~150 bp) with no nucleosome flanks called; runtime grows sharply on shallow data. Command: `macs3 hmmratac -i rep1.bam -f BAMPE -n rep1_hmm --outdir hmm/` (also writes a model .json and a cutoff-analysis table). The 30M-read minimum below and any runtime figures are practitioner guidelines, not tested here.

**Fix:** Verify fragment-size periodicity in QC first (atac-qc skill). If flat, fall back to MACS3 callpeak. HMMRATAC needs deep coverage (a practical minimum around 30M deduplicated nuclear reads, untested here).

### HOMER findPeaks -- Window-size sensitivity

**Trigger:** Default `-style dnase` uses 75 bp peaks; ATAC peaks are 250-500 bp typically.

**Mechanism:** HOMER's window-based caller does not auto-fit width to ATAC.

**Fix:** Use `-style factor -size 150` for narrow ATAC peaks, or skip HOMER for peak calling and use it only for downstream motif analysis on MACS peaks.

### Aligner choice -- chromap vs bwa-mem2 vs bowtie2 affects peak shape

**Trigger:** Switching aligners between datasets and expecting reproducible peaks.

**Mechanism:** chromap (Zhang 2021) applies its own ATAC-specific 4 bp / -5 bp Tn5 shift before fragment output; bwa-mem2 and bowtie2 do not. Downstream `--shift -75 --extsize 150` parameters are calibrated for unshifted bwa/bowtie BAMs; applying them to chromap output double-shifts the signal.

**Symptom:** Peaks called from chromap output are shifted by ~5-10 bp relative to bwa output at the same locus.

**Fix:** When using chromap, drop `--shift` and `--extsize` (chromap's pre-shift is sufficient) OR omit chromap's `--Tn5-shift` (the shift is opt-in, applied only when that flag or an ATAC preset is set) so it is not double-applied, then proceed with standard MACS parameters. Document the aligner version and any shift choices in methods. Within a project, pin the aligner.

### Single-sample (no replicate)

**Trigger:** Single biological sample without any replicate for IDR.

**Mechanism:** IDR requires two replicates by construction. For n=1, statistical confidence per peak comes from local background (Poisson p-value) but reproducibility cannot be assessed.

**Fix:** Call with `-q 0.05` (instead of the `-p 0.01` + IDR pattern) and report that reproducibility was not assessed; this is a single-sample setting, not ENCODE-compliant.

## ENCODE 3 vs ENCODE 4 Differences

| Feature | ENCODE 3 (legacy) | ENCODE 4 (current) |
|---------|-------------------|---------------------|
| Per-rep significance threshold | `-q 0.05` directly | `-p 0.01` (loose) + IDR |
| Pseudoreplicate IDR cutoff | Not formalized | same 0.05 threshold as true replicates; rescue and self-consistency ratios |
| TSS enrichment threshold | >= 6 (older) | >= 7 (hg38, GENCODE v29) |
| Mt fraction expectation | < 25% | < 20% (Omni-ATAC < 5%) |
| Blacklist | v1 | v2 (Amemiya 2019) |
| Default genome size | hardcoded `hs`/`mm` | encouraged: deepTools effectiveGenomeSize |

To reproduce a published ENCODE 3 dataset, pin the original pipeline and threshold exactly. ENCODE 4 results are not directly numerically comparable to ENCODE 3 even on the same input BAM.

## Super-Enhancers

Super-enhancer calling (ROSE; Whyte 2013) is defined on H3K27ac, MED1 or BRD4 ChIP-seq signal, not on ATAC-seq. ROSE-style stitching of ATAC peaks can be applied when no ChIP data exist, but it is not a standard ATAC-seq analysis, and its calls are not comparable to published super-enhancer sets. For the full workflow (stitching and TSS-exclusion parameters, signal ranking and inflection, marker choice, hg38 implementation caveats), use the `bio-chipseq-super-enhancers` Skill.

## ENCODE-Style ATAC-seq Pipeline (Reference Implementation)

An ENCODE-style pattern (`scripts/call_atac_peaks.sh` runs the whole pipeline). The ENCODE pipeline instead calls peaks on Tn5-shifted tagAlign (`-f BED`) with `--call-summits`, so peak sets are comparable but not identical. `$GSIZE` is the read-length-matched effective size (see above).

```bash
# Per-replicate peak calling (loose threshold)
macs3 callpeak \
    -t rep1.filt.dedup.bam \
    -f BAM -g $GSIZE \
    -n rep1 --outdir peaks/rep1/ \
    --nomodel --shift -75 --extsize 150 \
    --keep-dup all \
    -B --SPMR \
    -p 0.01

# Pooled (all replicates)
macs3 callpeak \
    -t rep1.filt.dedup.bam rep2.filt.dedup.bam \
    -f BAM -g $GSIZE -n pooled --outdir peaks/pooled/ \
    --nomodel --shift -75 --extsize 150 --keep-dup all -B --SPMR -p 0.01

# Pseudoreplicates: disjoint halves. -s picks reads by name hash (seed.fraction), -U writes the complement,
# so no read is in both halves. Two independent -s draws with different seeds would share ~50% of reads.
samtools view -b -s 1.5 -U rep1.psr2.bam -o rep1.psr1.bam rep1.filt.dedup.bam
# (call peaks on each pseudoreplicate the same way)
```

`macs2 callpeak` takes the same flags (see `usage-guide.md` for its import failure).

`--SPMR` writes signal as Signal Per Million Reads (normalized bedGraph). `-p 0.01` is intentionally loose; IDR will tighten to a reproducible set.

## IDR for Reproducible Peaks

**Goal:** Find peaks reproducible across biological replicates at controlled IDR.

**Approach:** Rank paired peak lists by p-value, fit IDR's two-component mixture (reproducible + noise), threshold at IDR <= 0.05 for every comparison below.

```bash
# Sort peaks by p-value (column 8) so IDR scores by significance
sort -k8,8nr rep1_peaks.narrowPeak > rep1.sorted.narrowPeak
sort -k8,8nr rep2_peaks.narrowPeak > rep2.sorted.narrowPeak

# True replicates -- threshold IDR <= 0.05
idr --samples rep1.sorted.narrowPeak rep2.sorted.narrowPeak \
    --input-file-type narrowPeak --rank p.value \
    --output-file true_reps.idr \
    --idr-threshold 0.05 --plot --log-output-file idr.log

# Pseudoreplicates -- same threshold (the script runs rep1, rep2 and pooled pseudoreplicate pairs)
idr --samples psr1_peaks.narrowPeak psr2_peaks.narrowPeak \
    --input-file-type narrowPeak --rank p.value \
    --output-file psr.idr --idr-threshold 0.05 --plot
```

**ENCODE reproducibility rules** (ENCODE-DCC/atac-seq-pipeline `encode_task_reproducibility.py`). Nt = peaks passing IDR on the true replicate pair; N1, N2 = on each replicate's own pseudoreplicate pair; Np = on the pooled pseudoreplicates. Rescue ratio = max(Np, Nt) / min(Np, Nt); self-consistency ratio = max(N1, N2) / min(N1, N2). Pass if both <= 2, borderline if one > 2, fail if both > 2. IDR score >= 540 corresponds to IDR <= 0.05 (score = int(-125 log2 IDR)).

**IDR fails when:** Ranking column choice matters. `--rank p.value` (column 8) is robust; `--rank signal.value` (column 7) breaks if MACS pile-up scaling differs between replicates.

## Decision Tree by Experimental Scenario

| Scenario | Recommended caller | Why |
|----------|-------------------|-----|
| Bulk ATAC, 2-3 reps, depth >= 25M | MACS3/MACS2 ENCODE-style pipeline + IDR | Reproducible, comparable to published peaksets |
| Bulk ATAC, 1 sample (no rep) | MACS3 callpeak with `-q 0.05`; do not run IDR | IDR is meaningless without reps; tighter q-value substitutes |
| Bulk ATAC, deep library, want NFR + flanking nuc structure | MACS3 hmmratac | HMM separates NFR from nucleosome flanks |
| Multi-replicate joint analysis where rep weighting is symmetric | Genrich `-j` ATAC mode | Joint p-value across reps; built-in chrM and blacklist |
| Cell type with broad accessible domains | MACS3 `--broad --broad-cutoff 0.1` for domains; narrow for differential | Domain-level inference vs site-level; super-enhancer calling is a separate ChIP-based workflow (see Super-Enhancers) |
| FFPE / degraded chromatin (flat fragment dist) | MACS3 callpeak with stringent `-q 0.01`; never HMMRATAC | HMM needs fragment periodicity |
| scATAC pseudobulk per cluster | MACS3 callpeak per cluster + iterative overlap | See atac-seq/single-cell-atac |
| Want fixed-width consensus peaks for differential | Call broadly, then re-center to summits +/- 250 bp | See atac-seq/consensus-peakset |
| Plant / non-model organism | MACS3 with `-g <effective_size>`; verify size empirically | Default `-g hs/mm` invalid; compute via khmer |

## Reconciliation: When Callers Disagree

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| MACS narrow peaks much fewer than Genrich | Different statistical models and peak merging (GM12878 chr1 slice: Genrich 14,328 joint peaks at `-q 0.05` vs MACS3 4,803 for one replicate; `-q 0.01` gave 22,936, so lowering `-q` did not reduce the count) | Do not tune `-q` for parity; compare by overlap with IDR peaks or keep the MACS-supported subset |
| HMMRATAC misses peaks MACS finds | Library too shallow OR fragment-size periodicity weak | Trust MACS; HMMRATAC is depth-sensitive |
| HMMRATAC calls peaks MACS misses | HMM is sensitive to mid-strength accessibility flanked by phased nucleosomes | Inspect; often genuine but unconfirmed by short-fragment signal |
| Same peak called by all but width 2x different | Broad mode vs narrow mode mismatch | Standardize: re-center to summit +/- 250 bp for differential |
| Per-rep MACS calls peak; pooled MACS does not | One rep dominates; lambda smoothes it out in pooled | Trust pooled + IDR over per-rep counts |

**Operational rule for high-confidence reporting:** Require a peak to pass IDR <= 0.05 on true replicates AND survive blacklist/greylist filtering AND have mean signalValue >= 5 across reps. Two callers from different families (MACS + Genrich) agreeing within 250 bp is acceptable evidence when IDR is unavailable.

## Blacklist and Greylist

```bash
# ENCODE blacklist (Amemiya 2019) -- always remove
wget https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg38-blacklist.v2.bed.gz
gunzip hg38-blacklist.v2.bed.gz
bedtools intersect -v -a peaks.narrowPeak -b hg38-blacklist.v2.bed > peaks.no_blacklist.narrowPeak

# Sample-specific greylist (input-derived high-signal regions; rarely available for ATAC)
# For ATAC, ENCODE recommends pooling all samples' top-percentile signal and removing
# regions exceeding 100x median coverage as a "soft greylist"
```

Blacklist is mandatory; greylist is optional and most useful when the same library prep produces consistent artifact regions across samples.

## NFR-Only Peak Calling

**Goal:** Call peaks using only sub-nucleosomal fragments (<100 bp) for sharper TF-binding-relevant accessibility.

**Approach:** Pre-filter BAM to short fragments, then call peaks with parameters scaled to the smaller fragment length.

```bash
samtools view -h sample.dedup.bam | \
    awk 'substr($0,1,1)=="@" || ($9 > 0 && $9 < 100) || ($9 < 0 && $9 > -100)' | \
    samtools view -b > nfr.bam
samtools index nfr.bam

macs3 callpeak -t nfr.bam -f BAM -g $GSIZE -n sample_nfr \
    --nomodel --shift -37 --extsize 75 \
    --keep-dup all -p 0.01
```

`--shift -37 --extsize 75` halves both parameters to match shorter fragments; this is a fragment-scaled convention for NFR-focused input, not a TOBIAS-specified setting (TOBIAS instead applies the +4/-5 Tn5 correction to the full BAM via ATACorrect).

## Output Files (narrowPeak)

| Column | Field | Notes |
|--------|-------|-------|
| 1-3 | chrom, start, end | 0-based, half-open |
| 4 | name | MACS auto-numbers |
| 5 | score | Min(int(-10*log10(qvalue)), 1000) |
| 6 | strand | `.` for ATAC |
| 7 | signalValue | Fold enrichment over local lambda |
| 8 | pValue | -log10 p |
| 9 | qValue | -log10 q (BH-FDR) |
| 10 | summit_offset | Peak summit relative to start |

Convert to bigWig for browsers: `sort -k1,1 -k2,2n sample_treat_pileup.bdg > sample.sorted.bdg && bedGraphToBigWig sample.sorted.bdg chrom.sizes sample.bw` (bedGraphToBigWig is multi-pass and cannot read from a pipe/stdin, so sort to a file first).

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| `--shift/--extsize ignored` warning | Used `-f BAMPE` with these flags | Switch to `-f BAM` or remove the flags |
| 0 peaks called | Forgot `--nomodel`; MACS tries to build a shifting model and fails | Add `--nomodel --shift -75 --extsize 150` |
| Peak count >> 500k | Did not deduplicate; or did not remove chrM; or `-q` too loose | Pre-filter (samtools view -F 1804 -q 30; samtools idxstats); use `-q 0.01` |
| `Sequence chrM not found` (Genrich) | Wrong chromosome name in `-e` flag (chrM vs MT) | Match BAM header naming convention |
| HMMRATAC out of memory | Deep library / large genome; the current tool is `macs3 hmmratac` (Python, not Java) | Increase available RAM and use a scratch `--outdir`; the `-Xmx`/`HMMRATAC.jar` heap flags apply only to the deprecated standalone Java HMMRATAC |
| Peaks shifted by 75 bp from expected positions | Forgot `--shift -75` (cuts at one end of read) | Add the shift; positions are now centered on Tn5 cut site |
| IDR returns 0 reproducible peaks | Sorted by wrong column; ranks are random | Sort each peakset by `-k8,8nr` (p-value descending) |

## References

- Buenrostro JD et al 2013 Nat Methods 10:1213 (ATAC-seq protocol)
- Corces MR et al 2017 Nat Methods 14:959 (Omni-ATAC protocol)
- Corces MR et al 2018 Science 362:eaav1898 (iterative-overlap fixed-width 501 bp consensus peaks)
- Tarbell ED & Liu T 2019 Nucleic Acids Res 47:e91 (HMMRATAC)
- Gaspar JM, Genrich: detecting sites of genomic enrichment (github.com/jsh58/Genrich; no published paper)
- Li Q et al 2011 Ann Appl Stat 5:1752 (IDR framework)
- Landt SG et al 2012 Genome Res 22:1813 (ENCODE/modENCODE peak calling guidelines, IDR Nself rule)
- Amemiya HM et al 2019 Sci Rep 9:9354 (ENCODE blacklist v2)
- ENCODE ATAC-seq Standards (encodeproject.org/atac-seq) -- canonical pipeline parameters
