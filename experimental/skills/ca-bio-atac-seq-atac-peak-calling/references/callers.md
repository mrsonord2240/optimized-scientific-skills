# Callers



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



```bash
# Genrich needs name-sorted BAMs (coordinate-sorted input aborts). BAMs are already deduplicated, so no -r.
samtools sort -n -o rep1.nsort.bam rep1.filt.dedup.bam
samtools sort -n -o rep2.nsort.bam rep2.filt.dedup.bam
Genrich -t rep1.nsort.bam,rep2.nsort.bam -o joint.narrowPeak -j -e chrM -E hg38-blacklist.v2.bed -q 0.05

# hmmratac writes <name>_accessible_regions.narrowPeak, a model .json and a cutoff-analysis table
macs3 hmmratac -i rep1.filt.dedup.bam -f BAMPE -n rep1_hmm --outdir hmm/
```
