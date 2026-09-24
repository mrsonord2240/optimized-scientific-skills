# Library Metrics Reference

Read this reference when evaluating insert-size distributions or GC bias. Use the main SKILL.md for scope, thresholds, validator exits, and cross-skill routing.

## Insert Size Distribution

**Goal:** Verify that the fragment length distribution matches the library preparation protocol.

**Approach:** Extract template_length from properly paired reads and compare the distribution to expected values for the library type.

### samtools stats

~~~bash
samtools stats input.bam > stats.txt
grep "^IS" stats.txt | cut -f2,3 > insert_sizes.txt
~~~

### Picard CollectInsertSizeMetrics

~~~bash
picard CollectInsertSizeMetrics \
    I=input.bam \
    O=insert_metrics.txt \
    H=insert_histogram.pdf
~~~

H= asks Picard to invoke Rscript for the histogram; omit H= if R is unavailable and retain O= for the metrics table.

### Expected Insert Sizes by Library

| Library | Mean insert | Distribution shape | Diagnostic |
|---------|-------------|-------------------|-----------|
| TruSeq DNA PCR-free WGS | 400-500 bp | Roughly Gaussian, tight/sharp | Sharpest, most symmetric peak (no PCR bias); bimodality = degraded sample |
| TruSeq DNA Nano (PCR) WGS | 300-400 bp | Gaussian, broadened by PCR amplification bias | |
| Twist / IDT exome capture | 250-350 bp | Gaussian | |
| TruSeq Stranded mRNA | 200-300 bp | Right-skewed (transcript distribution) | Long tail = poor size selection |
| Ribo-Zero rRNA-depleted | 250-400 bp | Right-skewed | |
| Smart-seq2 / Smart-seq3 | 200-700 bp | Broad | |
| 10x Chromium (3') | n/a -- not informative | n/a | |
| TruSeq ChIP | 200-400 bp | Sharp | |
| ATAC-seq (Buenrostro / Omni-ATAC) | Multimodal | Peaks at ~50, ~180, ~370 bp | **Missing multimodal pattern = bad library**; missing ~180 bp mononucleosome peak = Tn5 over- or under-titrated, or degraded DNA. Check with the samtools stats IS histogram above |
| Hi-C / Micro-C | Multimodal | Peak at ligation-junction size | |
| cfDNA / ctDNA | 160-180 bp | Multimodal; ~167 bp mononucleosomal + ~340 dinuc | Tumor-derived shorter (~145 bp); shape itself is a biomarker |
| FFPE | 100-250 bp | Right-skewed, broad | |
| aDNA | 30-80 bp | Sharp left-skewed | |
| ONT (native) | 1-30 kb | n/a | |
| PacBio HiFi | 10-25 kb | Sharp peak | |

### Python Insert Size Analysis

~~~python
import pysam
import numpy as np
import matplotlib.pyplot as plt

def get_insert_sizes(bam_file, max_reads=100000):
    sizes = []
    bam = pysam.AlignmentFile(bam_file, 'rb')
    for i, read in enumerate(bam.fetch(until_eof=True)):  # until_eof: no index needed
        if i >= max_reads:
            break
        if read.is_proper_pair and not read.is_secondary and read.template_length > 0:
            sizes.append(read.template_length)
    bam.close()
    return sizes

sizes = get_insert_sizes('sample.bam')
print(f'Median insert size: {np.median(sizes):.0f}')
print(f'Mean insert size: {np.mean(sizes):.0f}')
print(f'Std dev: {np.std(sizes):.0f}')

plt.hist(sizes, bins=100, range=(0, 1000))
plt.xlabel('Insert Size')
plt.ylabel('Count')
plt.savefig('insert_size_dist.pdf')
~~~

## GC Bias

GC content correlation with coverage.

### Picard CollectGcBiasMetrics

~~~bash
picard CollectGcBiasMetrics \
    I=input.bam \
    O=gc_bias_metrics.txt \
    CHART=gc_bias_chart.pdf \
    S=gc_summary.txt \
    R=reference.fa
~~~

CHART= likewise needs Rscript; omit it if only the O= / S= tables are needed.

The summary (S=) has AT_DROPOUT and GC_DROPOUT (percent of reads lost in the low-AT / low-GC windows); the detail file (O=) has NORMALIZED_COVERAGE per GC bin (1.0 = mean coverage). The GC bias bands in the Quality Thresholds table are applied to NORMALIZED_COVERAGE of well-populated bins.

### deepTools computeGCBias

~~~bash
computeGCBias \
    -b input.bam \
    --effectiveGenomeSize 2913022398 \
    -g hg38.2bit \
    -o gc_bias.txt \
    --biasPlot gc_bias.pdf
~~~

### Interpret GC Bias

| Issue | Symptom |
|-------|---------|
| Under-representation | Low GC coverage drops |
| Over-representation | High GC coverage elevated |
| PCR bias | Strong correlation |
