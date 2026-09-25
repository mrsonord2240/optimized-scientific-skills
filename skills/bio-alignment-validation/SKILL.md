---
name: bio-alignment-validation
description: Validate BAM integrity, reference-dictionary identity, alignment quality, contamination, and sample swaps with insert size, pairing, GC bias, strand balance, and post-alignment metrics. Use before variant calling or quantification.
tool_type: mixed
primary_tool: samtools
license: MIT
author: GPTomics
---

## Version Compatibility

Reference examples tested with: matplotlib 3.8+, numpy 1.26+, picard 3.1+, pysam 0.22+, samtools 1.19+
Checked on: samtools 1.24, pysam 0.24.1, Picard 3.5.0, VerifyBamID2 2.0.3, somalier 0.3.5, deepTools 4.0.0, RSeQC 5.0.4
Install: `conda install -c bioconda samtools picard` and `pip install pysam numpy matplotlib`. Install R when requesting Picard charts (`H=` / `CHART=`): Picard invokes `Rscript` for those plots, while omitting the chart argument still writes the metrics table. `picard <Tool> KEY=VALUE` is the bioconda wrapper for `java -jar picard.jar <Tool> KEY=VALUE`; Picard 3.x still accepts the `KEY=VALUE` form.

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Alignment Validation

Post-alignment quality control to verify alignment quality and identify issues.

**"Check alignment quality"** -> Compute post-alignment QC metrics (mapping rate, pairing, insert size, strand balance) to identify issues before downstream analysis.
- CLI: `samtools flagstat`, `samtools stats`, Picard `CollectAlignmentSummaryMetrics`
- Python: `pysam.AlignmentFile` iteration with metric calculations

## What Each Validation Covers

## Reference Files

- `references/library-metrics.md` — Read when selecting or interpreting insert-size and GC-bias metrics, including the Picard chart/R requirement and the runnable Python insert-size example.

| Concern | Tools | What it catches |
|---------|-------|-----------------|
| **File integrity** | `samtools quickcheck`, `samtools view -c`, `picard ValidateSamFile` | See the detection table below: each tool sees a different slice of the damage |
| **Sequence dictionary identity** | `samtools dict` + SN/LN/M5 comparison | BAM aligned to wrong reference flavor / different decoy / chr vs no-chr |
| **QC metrics** | `samtools stats`, `flagstat`, `mosdepth`, Picard `CollectMultipleMetrics` / `CollectHsMetrics` / `CollectWgsMetrics`, `examples/validate_alignment.*`; `references/library-metrics.md` | Are the data biologically reasonable for the assay? Read the reference for insert size and GC bias. |
| **Contamination / sample swap** | `verifybamid2`, `somalier`, Picard `CrosscheckFingerprints` | Cross-sample contamination, tumor-normal swap, mislabeled sample |

A file can pass `quickcheck` and still be malformed in ways that crash GATK three hours into HaplotypeCaller. Conversely, a QC-poor BAM can be structurally valid.

### File Integrity
```bash
# Fast: header + EOF block check (misses mid-file truncation, invalid CIGAR)
samtools quickcheck -v in.bam || echo "QUICKCHECK FAILED"
samtools quickcheck -v *.bam > bad_bams.fofn   # one fail-line per bad file

# Slow but thorough: structural validation (R= enables the NM check; use samtools calmd to check MD)
picard ValidateSamFile I=in.bam MODE=SUMMARY R=ref.fa
```

What each check caught in 21 planted defects (synthetic BAMs derived from a real human PE BAM; samtools 1.24, Picard 3.5.0), and no flags on the 2 valid synthetic controls:

| Check | Caught | Misses |
|-------|--------|--------|
| `samtools quickcheck` | 3/21: missing EOF block, truncated file, header without `@SQ` | everything inside the file: bit flips, dropped blocks, bad records |
| `samtools view -c` (decodes every record) | 3/21: truncated file, bit-flipped block, CIGAR/SEQ length mismatch | dropped blocks and every semantic defect |
| `picard ValidateSamFile` | 18/21: mate fields (`MATE_NOT_FOUND`, `MISMATCH_MATE_ALIGNMENT_START`, `MISMATCH_FLAG_MATE_*`), unsorted (`RECORD_OUT_OF_ORDER`), read group missing from header, MAPQ on unmapped reads, CIGAR off the contig end, CIGAR/SEQ length, bit flips, a dropped mid-file block (as orphaned mates) | 0x40 and 0x80 both set, TLEN inconsistent with mate positions, header-only file with 0 records |
| CI one-liner below | 6/21: the quickcheck and decode catches plus 0-record files | all record-level defects |

No tool here sees a whole block removed from the middle unless mates are orphaned; for transfers compare `md5sum` or the read count with the source.

Picard noise on valid files (observed): `MATE_NOT_FOUND` on region-extracted BAMs (mates fall outside the slice), `MISSING_TAG_NM` warnings on BAMs without NM tags (STAR output), `RECORD_OUT_OF_ORDER` on a samtools name-sorted BAM (Picard likely orders query names differently from samtools). Long-read BAMs can also report missing `@RG PL` / other required header tags and NM-convention errors. Read the SUMMARY before trusting a non-zero count.

`IGNORE=<TYPE>` removes a check. On the planted files `IGNORE=INVALID_MAPPING_QUALITY IGNORE=MISMATCH_FLAG_MATE_NEG_STRAND` turned 40, 2 and 2820 genuine errors into "No errors found". Ignore a type only after tracing it to a known harmless source.

CI-safe one-liner (`MIN_READS` = smallest read count you accept; default 1):
```bash
MIN_READS=${MIN_READS:-1}
test -s in.bam \
  && samtools quickcheck -v in.bam \
  && [ "$(samtools view -c -F 2304 in.bam)" -ge "$MIN_READS" ] \
  || { echo "BAM failed integrity"; exit 1; }
```

### Sequence Dictionary Cross-Validation (SN / LN / M5)
```bash
# Compare contig names and lengths, and M5 wherever the BAM header carries it. Exit 1 on any difference.
awk -F'\t' '
  function parse(   i,k,v) { sn=ln=m5=""; for (i=2;i<=NF;i++) { k=substr($i,1,2); v=substr($i,4); if (k=="SN") sn=v; else if (k=="LN") ln=v; else if (k=="M5") m5=v } }
  FNR==NR { if ($1=="@SQ") { parse(); rlen[sn]=ln; rm5[sn]=m5 } ; next }
  $1=="@SQ" { parse(); if (m5!="") nm5++
      if (!(sn in rlen))      { print "NOT IN REFERENCE:", sn; bad=1 }
      else if (rlen[sn]!=ln)  { print "LENGTH DIFFERS:", sn, ln, "vs", rlen[sn]; bad=1 }
      else if (m5!="" && m5!=rm5[sn]) { print "M5 DIFFERS:", sn; bad=1 } }
  END { if (!nm5) print "no M5 in BAM header: only names and lengths were compared"; exit bad }
' <(samtools dict ref.fa) <(samtools view -H in.bam)
```

Names and lengths (SN, LN) catch renaming (UCSC `chr1` vs Ensembl `1`), a different contig set and different lengths. M5 catches a different sequence under the same name and length, and a contig-to-sequence mix-up. Most aligners (bwa, minimap2, STAR) write no M5, so a BAM without it can only be checked on SN/LN; the snippet says so instead of failing. Concrete failure modes: GRCh38 vs GRCh38.p13 vs GRCh38_no_alt (alt contigs differ); soft-masked vs unmasked (M5 matches -- case is normalized to uppercase before hashing, so lowercase soft-masking is invisible to M5). Hard-masking is different: it replaces bases with `N`, changing the sequence, so its M5 does NOT match the unmasked reference. Where M5 is present it is the only definitive identity check.

### Contamination and Sample Swap

No alignment QC is complete without these in production. `verifybamid2` and `somalier` need whole-genome or exome-scale data covering their marker panels (a region slice leaves too few markers: VerifyBamID2 2.0.3 reported "Insufficient Available markers" on a 100 kb slice).
```bash
# Cross-sample contamination; the prefix of the released SVD files includes ".dat"
verifybamid2 --SVDPrefix /resources/1000g.phase3.10k.b38.vcf.gz.dat \
    --Reference ref.fa --BamFile sample.bam --Output sample.contam
# FREEMIX > 0.03 is the commonly used contamination-concern threshold; values escalate from there.

# Relatedness, sex check, sample swap detection
somalier extract -d extracted/ -s /resources/sites.hg38.vcf.gz \
    -f ref.fa sample.bam
somalier relate --infer extracted/*.somalier

# Tumor/normal pairing verification: compare whole files so duplicate RG IDs/PUs cannot collapse
# two samples into one group. Alternatively give every BAM a distinct RG ID and PU.
# Give paired tumor/normal BAMs the same SM (patient id) or add EXPECT_ALL_GROUPS_TO_MATCH=true.
picard CrosscheckFingerprints I=tumor.bam I=normal.bam \
    HAPLOTYPE_MAP=Homo_sapiens_assembly38.haplotype_database.txt \
    CROSSCHECK_BY=FILE LOD_THRESHOLD=-5 OUTPUT=crosscheck.metrics
# With LOD_THRESHOLD=-5: LOD > 5 = same individual; < -5 = different; in between = ambiguous
```

Checked here: `CrosscheckFingerprints` and `somalier extract` + `relate --infer` ran end-to-end on a 40-site chr22 haplotype map / sites VCF built for the test slice (same reads: LOD +12.6, somalier relatedness 1.0; genotype-flipped reads: LOD -63, exit 1 for a same-SM pair, somalier relatedness -140, concordance 0). `verifybamid2` ran to the marker check only, so no FREEMIX value was produced here; the `HAPLOTYPE_MAP` above is the GATK-bundle file and the somalier/VerifyBamID2 panels are downloads, none shipped with this Skill.

## Proper Pairing Rate

Percentage of mapped paired reads that are correctly paired. Every count below is over primary records (`-F 2304`; mapped: `-F 2308`), the same denominators as `flagstat` "primary".

### samtools flagstat

```bash
samtools flagstat input.bam

samtools flagstat input.bam | grep "properly paired"
```

### Calculate Pairing Rate

```bash
proper=$(samtools view -c -f 2 -F 2308 input.bam)
paired=$(samtools view -c -f 1 -F 2308 input.bam)
awk -v p="$proper" -v n="$paired" 'BEGIN{if (n) printf "Proper pairing rate: %.2f%%\n", 100*p/n; else print "no paired reads"}'
```

Multiply before dividing: `bc` with `scale=2` on `$proper / $paired * 100` truncates the quotient first (99.96% prints as 99.00%).

## Strand Balance

The metric is the **forward fraction F/(F+R)** over mapped primary reads. A balanced 0.48-0.52 applies to WGS / WES / generic DNA-seq on autosomes (do not use the F/R quotient, which is ~1.0 when balanced). Expected to deviate for: stranded RNA-seq (deliberately strand-asymmetric -- verify with RSeQC `infer_experiment.py -r genes.bed12 -i rna.bam`, needs a BED12 gene model: ~0.5/0.5 = unstranded, one fraction near 1 = stranded; RSeQC 5.0.4), bisulfite (CT vs GA), small-RNA / strand-specific RNA-seq, and chrY/chrM regions. Per-chromosome strand imbalance >5% on autosomes is a field-convention rule of thumb (no single primary citation) — it picks up aligner artifacts; on chrX/chrY it suggests sex-mismatch.

### Calculate Strand Fraction

```bash
forward=$(samtools view -c -F 2324 input.bam)          # mapped, primary, forward
reverse=$(samtools view -c -f 16 -F 2308 input.bam)    # mapped, primary, reverse
echo "Forward: $forward"
echo "Reverse: $reverse"
awk -v f="$forward" -v r="$reverse" 'BEGIN{if (f+r) printf "Forward fraction: %.3f\n", f/(f+r)}'
```

### Check Strand Bias per Chromosome

Needs an index (`samtools index input.bam`).


```bash
samtools idxstats input.bam | awk '$1!="*" && $3>0 {print $1}' | head -25 | while read -r chr; do
    fwd=$(samtools view -c -F 2324 input.bam "$chr")
    rev=$(samtools view -c -f 16 -F 2308 input.bam "$chr")
    awk -v c="$chr" -v f="$fwd" -v r="$rev" 'BEGIN{if (f+r) printf "%s: F=%d R=%d forward fraction=%.3f\n", c, f, r, f/(f+r)}'
done
```

## Mapping Quality Distribution

Mapped primary reads only (`-F 2308`): unmapped records carry MAPQ 0 (or a stray value) and would distort both the histogram and the mean.

### Extract MAPQ Distribution

```bash
samtools view -F 2308 input.bam | cut -f5 | sort -n | uniq -c | sort -k2 -n
```

### Calculate Mean MAPQ

```bash
samtools view -F 2308 input.bam | awk '{sum+=$5; count++} END {if (count) print "Mean MAPQ:", sum/count}'
```

### MAPQ Distribution Is Bimodal and Aligner-Specific

Mean MAPQ is misleading; distributions are bimodal (0 and aligner-max). For aligner-specific scales and "unique mapping" sentinels, see bio-sam-bam-basics. The fraction of primary mapped reads at MAPQ >= 30 (`samtools view -c -F 2308 -q 30` over `-c -F 2308`) is a more informative summary than the mean.

## Chromosome Coverage Balance

### Calculate Per-Chromosome Coverage

Mapped reads per bp (a density, not depth):
```bash
samtools idxstats input.bam | awk '$2>0 {printf "%s\t%.4f\n", $1, $3/$2}' | head -25
```

### Check for Aneuploidy / Sex Chromosome Imbalance

Median-normalized per-autosome density (1.0 = expected diploid; 0.5 = monosomy/sex; 1.5 = trisomy). Plain POSIX awk (no gawk `asort`), plain `chr1`-`chr22` / `1`-`22` only (alt, decoy, HLA and unplaced contigs would distort the median):
```bash
samtools idxstats in.bam | awk '$2>0 && $3>0 && $1 ~ /^(chr)?[0-9]+$/ {print $1, $3/$2}' \
  | sort -k2,2g \
  | awk '{c[NR]=$1; v[NR]=$2} END {
      if (NR==0) { print "no autosomes with reads" > "/dev/stderr"; exit 1 }
      med = v[int(NR/2)+1]
      for (i=1; i<=NR; i++) printf "%s\t%.3f\n", c[i], v[i]/med }'
```

For full ancestry / contamination / relatedness checking, use `verifybamid2` or `somalier` -- they account for population AFs, not just per-contig depth.

## Mismatch Rate

### Picard CollectAlignmentSummaryMetrics

```bash
picard CollectAlignmentSummaryMetrics \
    I=input.bam \
    R=reference.fa \
    O=alignment_summary.txt
```

### Key Metrics

Picard reports these as fractions (0.9996 = 99.96%); the Good values are percentages, so multiply by 100.

| Metric | Description | Good Value |
|--------|-------------|------------|
| PCT_PF_READS_ALIGNED | Mapped fraction | > 95% (0.95) |
| PF_MISMATCH_RATE | Mismatches | < 1% (0.01) |
| PF_INDEL_RATE | Indels | < 0.1% (0.001) |
| STRAND_BALANCE | Strand ratio | ~0.5 |

## Comprehensive Validation Script

**Goal:** Run all key alignment QC checks in a single pass and generate a summary report.

**Approach:** `examples/validate_alignment.sh` runs `quickcheck`, then `flagstat`, `stats`, `idxstats` and primary-read counts, grades mapping rate, proper pairing, strand balance and mean MAPQ against the Quality Thresholds table, and prints the verdict.

```bash
bash examples/validate_alignment.sh sample.bam | tee qc_report.txt   # tee saves the report; exit status is the script's
```

Exit status: 0 = every metric PASS or WARN, 1 = a metric in the FAIL band, 2 = file missing, empty, unreadable or without primary records. The pipe to `tee` hides the script's status; use `set -o pipefail` or `${PIPESTATUS[0]}` in pipelines.

## Python Validation Module

`examples/validate_alignment.py` computes the same metrics with pysam and the same exit codes (0 / 1 / 2):
```bash
python examples/validate_alignment.py sample.bam        # whole file, no index needed
python examples/validate_alignment.py sample.bam -n 100000
```

It reads with `fetch(until_eof=True)`: unplaced unmapped reads (RNAME `*`) are counted, whereas `fetch()` and `get_index_statistics()` never see them (a BAM with 30% unplaced unmapped reads then reads as 100% mapped). `-n N` prints a head-of-file bias warning: it samples the first N records only, and in a coordinate-sorted BAM the unplaced tail is not reached, so the mapping rate is overestimated. Both validators require `@SQ` and do not grade pairing or strand balance below 100 primary records; they print that limit instead of interpreting a one-read fraction. For unbiased per-chromosome statistics use `samtools view -s 42.01 input.bam` (the `INT.FRAC` form uses `INT` as the seed; specify an explicit nonzero seed so the subsample is documented and consistent across paired runs).

## Quality Thresholds Summary

Germline short-read DNA defaults; the example validators grade mapping rate, proper pairing, strand balance and mean MAPQ with these bands. Per-assay values (RNA-seq, ATAC, ChIP, long-read, aDNA, ...) are in bio-bam-statistics.

| Metric | Good | Warning | Fail |
|--------|------|---------|------|
| Mapping rate | > 95% | 90-95% | < 90% |
| Proper pairing | > 90% | 80-90% | < 80% |
| Singletons (`flagstat`) | < 5% | 5-10% | > 10% |
| Duplicate rate (assay-specific) | see bio-bam-statistics decision table | -- | -- |
| Strand balance (forward fraction) | 0.48-0.52 | 0.45-0.55 | Outside |
| Mean MAPQ | > 40 | 30-40 | < 30 |
| GC bias (`NORMALIZED_COVERAGE`) | < 1.2x | 1.2-1.5x | > 1.5x (field-convention bands; Picard CollectGcBiasMetrics does not prescribe specific cutoffs) |

## When a Metric Fails

| Symptom | Check |
|---------|-------|
| Low mapping rate (< 90%) | Correct reference genome version; sample species matches reference; contamination with other species; read quality/adapters before alignment |
| Poor proper pairing (< 80%) | Insert size vs aligner expectation; chimeric reads from structural variants; library protocol; re-align with adjusted parameters |
| Abnormal insert size | Bimodal = mixed libraries; very wide = degraded DNA; adapter contamination; compare with the library table |
| Strand imbalance | Capture bias in targeted sequencing; PCR artifacts; look at specific chromosomes for local bias |
| GC bias | Strong correlation = PCR amplification issues, library complexity or low-input DNA; may need GC correction before analysis |
| Low MAPQ | Repetitive / duplicated regions with multi-mappers; longer reads or paired-end help; filter low MAPQ for variant calling |

## Related Skills

- bio-bam-statistics - Per-assay metric thresholds, depth/coverage tools, mosdepth
- bio-alignment-filtering - Aligner-specific MAPQ thresholds (canonical home)
- bio-duplicate-handling - Library-aware dedup decisions
- bio-sam-bam-basics - MAPQ-by-aligner table
- bio-chipseq-qc - ChIP-specific QC (FRiP, NSC, RSC)
