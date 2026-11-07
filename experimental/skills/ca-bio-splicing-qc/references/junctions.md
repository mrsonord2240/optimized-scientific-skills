## Junction Saturation

**Goal:** Determine whether sequencing depth is sufficient for comprehensive splicing detection.

**Approach:** RSeQC subsamples the splice events (5%, 10%, ... 100%) and counts junctions at each step. The numbers are only inside `*.junctionSaturation_plot.r` (vectors `x` percent, `y` known, `z` all, `w` novel); the helper parses them.

```bash
junction_saturation.py -i sample.bam -r genes.bed12 -o sample_junc_sat --skip-plot   # finer curve: lower -s and keep -l equal to -s (-l 2 -u 100 -s 2); -l != -s mislabels the percentages
python scripts/splicing_qc.py saturation sample.bam genes.bed12 sample_junc_sat
```

**Plateau rule:** if the **known**-junction curve grows by <2% from 80% to 100% of reads, it has plateaued; still rising means more sequencing would yield more known junctions. RSeQC shuffles the events without a seed, so the curve is stochastic (three repeats on one planted BAM differed by up to 16 known junctions at the same step; growth from 80% to 100% ranged 4.0-5.3%): compare libraries by growth percentage, not by exact counts.

A curve that is flat from the first steps is a *saturated* library, not an uninformative one; a curve of all zeros means no spliced reads or a BED/BAM contig mismatch.

## Novel-vs-Known Junction Ratio

**Goal:** Detect annotation/mapping issues or biologically interesting cryptic splicing.

```bash
junction_annotation.py -i sample.bam -r genes.bed12 -o sample_junc_annot --skip-plot
```

```python
import pandas as pd

# RSeQC .junction.xls: chrom, intron_st(0-based), intron_end(1-based), read_count, annotation
junc = pd.read_csv('sample_junc_annot.junction.xls', sep='\t')
junc['annotation'] = junc['annotation'].str.strip()    # RSeQC writes ' annotated' with a leading space
by_class = junc.groupby('annotation')['read_count'].sum()
known = by_class.get('annotated', 0) / by_class.sum()
novel = (by_class.get('partial_novel', 0) + by_class.get('complete_novel', 0)) / by_class.sum()
assert abs(known + novel - 1) < 1e-9, by_class.index.tolist()
print(f'known: {known:.1%}, novel: {novel:.1%}')
```

Without the `str.strip()` the class names never match and the snippet prints `known: 0.0%, novel: 0.0%` on every library.

**Definitions.** RSeQC calls a junction `annotated` when its donor AND acceptor are each in the gene model, so an unannotated skipping junction between two annotated exons counts as known; `partial_novel` has one known end, `complete_novel` none. The snippet above is **read-weighted**, which is what the table below uses. RSeQC stdout prints both **Splicing Events** (read-weighted) and **Splicing Junctions** (distinct-junction) blocks; do not apply the read-weighted table to the latter. On one planted library they were 93.4% and 34.6% known, respectively.

| Known fraction (reads) | Status | Interpretation |
|------------------------|--------|----------------|
| >=80% | Healthy | Comprehensive annotation, good alignment |
| 60-80% | Acceptable | Check annotation completeness or organism |
| <60% | Suspect or interesting | Mapping artifacts, contamination, OR biology |

If novel% >40%, first summarize support by annotation class; weak novel support or anchors suggest mapping noise, while strongly anchored and repeatedly supported novel sites warrant an annotation and biology check. This does not prove either cause:

```bash
# Uses the matching BAM and junction_annotation.py table. `reads` is unique,
# fragment-level support with >=8-nt anchors; zero-support rows are retained.
python scripts/splicing_qc.py novel-qc sample.bam sample_junc_annot.junction.xls --min-overhang 8
```

Then compare the candidate sites against the intended annotation release and inspect alignments. **High novel-junction rate may be biology, not artifact:**
- **TDP-43 loss** (ALS/FTD post-mortem brain): cryptic exon de-repression in UNC13A, STMN2, ATG4B (Brown 2022 *Nature*; Klim 2019 *Nat Neurosci*)
- **SF3B1-mutant** cancer (MDS, CLL, uveal melanoma): cryptic 3'ss ~10-30nt upstream of canonical (Darman 2015 *Cell Rep*)
- **Non-model organism**: GENCODE-grade annotation unavailable; novel junctions reflect annotation gaps

## Junction Read Overhang and Coverage

**Goal:** Per-junction read support and anchor lengths, to find weakly supported junctions.

```bash
python scripts/splicing_qc.py junctions sample.bam --min-overhang 8
```

`junction_stats()` in `scripts/splicing_qc.py` (works on an unindexed BAM) reports, per junction `(contig, intron_start_0based, intron_end)`: `reads` (fragments passing the overhang filter), `reads_all`, and `min_overhang`. Rules it applies, each checked on planted CIGARs:
- **Overhang** is the aligned length (M/=/X) of the block immediately left and right of each N, not the total matched on that side: `30M1000N4M800N66M` has overhang 4 at both junctions (summing gives 30).
- Secondary and supplementary records are skipped; NH>1 (or MAPQ < 30 without an NH tag) is skipped; the two mates of a pair count once per junction. On real chrX pass-2 data this equals STAR's unique-read `SJ.out.tab` count on all 2,765 shared junctions with `--min-overhang 0` (371 junctions with >=10 reads, same as STAR), while per-record counting agreed on only 2,182.
- `=`/`X` CIGAR operations are handled (an `=`-blind version put a junction 50 nt off).

Junction reads with overhang <8 nt are common false positives, especially for novel sites; most callers default to >=8 nt. For very deep BAMs use STAR `SJ.out.tab` or `regtools junctions extract` instead (the helper keeps read names in memory).
