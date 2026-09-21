---
name: bio-sashimi-plots
description: Creates sashimi-style plots showing RNA-seq read coverage and splice junction counts using ggsashimi (general-purpose, condition-grouped overlays), rmats2sashimiplot (rMATS-output-aware), MAJIQ-VOILA (LSV posteriors, interactive viewer; licence-gated), leafviz (leafcutter clusters Shiny), Jutils (tool-agnostic heatmaps and sashimi for rMATS/leafcutter/MntJULiP/MAJIQ output), or pyGenomeTracks (multi-track publication figures). Tool choice depends on the upstream differential-splicing tool's output format and the publication vs interactive use case. Use when visualizing specific splicing events, validating differential splicing calls, or producing publication-quality figures.
tool_type: python
primary_tool: ggsashimi
license: MIT
author: GPTomics
---

## Version Compatibility

Reference examples tested with (2026-09-20): ggsashimi 1.1.5, rmats2sashimiplot 4.0.0, leafcutter/leafviz 0.2.9, pyGenomeTracks 3.9, Jutils 1.5, pysam 0.24, pandas 2.2+, R 4.2.3 with **ggplot2 3.4.4**. MAJIQ/VOILA is licence-gated and was not installed: its commands below follow MAJIQ's public docs and were not run.

**ggsashimi needs ggplot2 < 3.5.** With ggplot2 3.5.2 and 4.0.3 the gene-model track and the per-group panels are shifted against each other and x tick labels are clipped (arc counts stay correct); with 3.4.4 everything aligns. **Always look at the rendered figure** for exon edges lining up with arcs before reporting it.

Install (ggsashimi and Jutils are on neither conda nor PyPI; PyPI `jutils` is an unrelated package):

```bash
conda install -c conda-forge -c bioconda rmats2sashimiplot pygenometracks pysam bedtools regtools samtools
conda install -c conda-forge r-base=4.2 r-ggplot2=3.4.4 r-data.table r-gridextra r-gtable seaborn scikit-learn
git clone https://github.com/guigolab/ggsashimi     # ggsashimi.py: one script, needs pysam and R on PATH
git clone https://github.com/splicebox/Jutils       # jutils.py
```

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Sashimi Plot Visualization

Visualize RNA-seq coverage tracks with splice junction arcs labeled by read count. Sashimi plots originated with MISO (Katz 2010 *Nat Methods*); modern tools differ in input handling, group aggregation logic, and customization. Tool choice is not interchangeable — some tools work only with specific upstream output formats.

## Tool Selection Matrix

| Tool | Best for | Input | Strengths | Fails when |
|------|----------|-------|-----------|------------|
| ggsashimi | Publication-quality grouped overlays from any BAM | BAMs + region | `--overlay` aggregates samples within a group; clean PDFs | No native rMATS/MAJIQ integration; need to extract coords manually; ggplot2 >= 3.5 misaligns panels |
| rmats2sashimiplot | One-line plot from rMATS output | rMATS event file + BAMs | No manual coord extraction | rMATS-specific; doesn't handle leafcutter or MAJIQ |
| MAJIQ-VOILA | Interactive LSV browsing with posterior PSI distributions | MAJIQ build + psi/deltapsi | Splice-graph topology; LSV-aware; posterior violins | Static figures; licence-gated download, not run in testing |
| leafviz | Cluster-level interactive browsing | leafcutter differential output | Filter table + sashimi-like plots | leafcutter-specific |
| Jutils | Unified output across rMATS, leafcutter, MntJULiP, MAJIQ | Tool-specific differential output | Heatmaps, Venn, sashimi tool-agnostically | Output less polished than ggsashimi |
| pyGenomeTracks | Multi-track publication figures (RNA-seq + ChIP/ATAC) | bedGraph or BigWig + BED + GTF | Combine RNA with chromatin tracks | Not splicing-specific; configure tracks manually |
| IGV (interactive) | Quick ad-hoc inspection | BAM + region | Scrollable, instant | Not for publication figures |
| MISO sashimi | Historical | MISO output | Original sashimi format | MISO unmaintained; no longer recommended |

## Decision Tree by Goal

| Goal | Recommended tool |
|------|-------------------|
| Validate a specific rMATS hit | rmats2sashimiplot (one-line) or ggsashimi (custom) |
| Validate a leafcutter cluster | leafviz (interactive) or ggsashimi with cluster coordinates |
| Validate a MAJIQ LSV (complex topology) | MAJIQ-VOILA (only tool that shows full LSV graph) |
| Publication-quality two-condition comparison | ggsashimi `-O 3 -A mean_j` for grouped overlay |
| Multi-track figure (RNA-seq + H3K4me3 + ATAC) | pyGenomeTracks |
| Quick ad-hoc browsing during development | IGV sashimi |
| Tool-agnostic batch heatmap of significant events | Jutils |
| Interactive cohort-level filtering of leafcutter results | leafviz Shiny |

## ggsashimi for Publication Overlays

**Goal:** Generate publication-quality sashimi plot for a region with samples grouped by condition and per-sample tracks aggregated.

**Approach:** Define samples + groups in a TSV (no header), a palette file, then call ggsashimi with coordinates, GTF, and visual flags. ggsashimi exits 0 when it drops a missing BAM, draws an empty region or hits an R error, so check the inputs before and the figure after.

```python
import subprocess
import pandas as pd
from pathlib import Path

# ggsashimi input: col1 = sample id, col2 = BAM path, col3 = group (used by -O overlay and -C colour)
groups = pd.DataFrame({
    'sample_id': ['ctrl1', 'ctrl2', 'ctrl3', 'trt1', 'trt2', 'trt3'],
    'bam': ['ctrl1.bam', 'ctrl2.bam', 'ctrl3.bam', 'trt1.bam', 'trt2.bam', 'trt3.bam'],
    'group': ['Control', 'Control', 'Control', 'Treatment', 'Treatment', 'Treatment']
})
groups.to_csv('sashimi_groups.tsv', sep='\t', index=False, header=False)
Path('palette.txt').write_text('#1f77b4\n#ff7f0e\n')  # one colour per group, in order of first appearance

missing = [b for b in groups['bam'] if not Path(b).is_file()]
assert not missing, f'ggsashimi would drop these BAMs silently: {missing}'

subprocess.run([
    'ggsashimi.py',
    '-b', 'sashimi_groups.tsv',
    '-c', 'chr17:43094000-43125000',   # contig spelled as in the BAM header
    '-o', 'BRCA1_sashimi',
    '--alpha', '0.25',
    '--height', '3',
    '--width', '10',
    '--shrink',
    '--fix-y-scale',
    '--ann-height', '4',
    '-g', 'gencode_v45.gtf',
    '--base-size', '14',
    '-O', '3', '-C', '3', '-P', 'palette.txt',
    '-A', 'mean_j',
    '-F', 'pdf'
], check=True)
assert Path('BRCA1_sashimi.pdf').is_file() and Path('BRCA1_sashimi.pdf').stat().st_size > 0, 'no figure written (R error above?)'
```

`examples/plot_sashimi.py` wraps this as `plot_sashimi()` and adds contig mapping (chrX vs X), an empty-region check and the `--shrink` guard below. Checked: ggsashimi's junction labels equal an independent pysam count (planted 3v3 set, real ENCODE 12-BAM locus, real chrX BAMs).

Key ggsashimi flags (Garrido-Martin 2018 *PLoS Comput Biol*):
- `-O 3`: column 3 of the TSV is the overlay level; samples of a group are drawn in one track. Required for `-A`
- `-C 3 -P palette.txt`: colour by column 3 using the palette file (R colour names or hex, one per line). Without `-C` everything is grey; with `-C` and no `-P` the colours are R defaults (red/green), not blue/orange
- `-A mean_j`: the arc label is the rounded (half to even: 6.5 shows 6) plain mean of the raw junction counts of the group's samples that have the junction; every sample's coverage stays overlaid. `mean` and `median` also aggregate the coverage; `median_j` is the median analogue. There is no depth normalization
- `-M`: minimum reads for a junction to be drawn, **default 1**, inclusive, applied **per sample before `-A`**: samples below `-M` drop out of the mean, biasing labels upward (planted 10.33 shows 11 at `-M 10`; 10 at `-M 1`). Keep `-M 1` for exact labels; raise it (5-10, 20+ for crowded plots) only to declutter, and never so high that no junction passes
- `--shrink`: rescale long introns for compact display. Crashes with `RuntimeError: generator raised StopIteration` when no junction passes `-M` (with `-s`, when one strand has none): lower `-M` or drop `--shrink`
- `--fix-y-scale`: identical y-axis across groups (essential for visual comparison)
- `--alpha 0.25`: transparency of per-sample coverage in overlay mode
- `--height`/`--width`/`--ann-height`/`--base-size`: sizes in inches / font size (e.g. `--width 12 --height 4 --base-size 16`)
- `-F pdf|svg|png|jpeg|tiff` (`-R` for raster PPI); pick explicitly: PDF for publication, PNG for slides, SVG for editing
- `-g`: GTF with exons; ggsashimi has no feature filter, so pre-filter the GTF (`awk '$0 ~ /protein_coding/'`) to restrict it
- Colour convention: control = blue (`#1f77b4`), treatment = orange (`#ff7f0e`); document it, use ColorBrewer for >2 groups

## Batch Plotting from rMATS Hits

**Goal:** Auto-generate sashimi plots for all significant rMATS differential events.

**Approach:** Parse SE.MATS.JC.txt, expand coordinates to flanking exons + 500nt context, map the contig name onto the BAM header, iterate ggsashimi and check every figure. rMATS writes `chrX`; an Ensembl-style BAM calls it `X` and ggsashimi dies with `ValueError: invalid contig`. Events near a contig start give a start < 1.

```python
import re
import subprocess
import pandas as pd
import pysam
from pathlib import Path

contigs = set(pysam.AlignmentFile(groups['bam'][0]).references)  # groups = the TSV above

def bam_contig(name):
    for cand in (name, name.removeprefix('chr'), 'chr' + name.removeprefix('chr')):
        if cand in contigs:
            return cand
    raise ValueError(f'contig {name} not in the BAM header')

diff = pd.read_csv('rmats_output/SE.MATS.JC.txt', sep='\t')
sig = diff[(diff['FDR'] < 0.05) & (diff['IncLevelDifference'].abs() > 0.10)]

Path('sashimi_plots').mkdir(exist_ok=True)
failed = []
for _, ev in sig.head(25).iterrows():
    region = f'{bam_contig(ev["chr"])}:{max(1, ev["upstreamES"] - 500)}-{ev["downstreamEE"] + 500}'
    safe_name = re.sub(r'[^A-Za-z0-9._-]', '_', f'{ev["geneSymbol"]}_{ev["chr"]}_{ev["upstreamES"]}_{ev["ID"]}')
    out = Path(f'sashimi_plots/{safe_name}.pdf')
    subprocess.run([
        'ggsashimi.py', '-b', 'sashimi_groups.tsv', '-c', region,
        '-o', str(out.with_suffix('')), '-M', '1', '--shrink', '--fix-y-scale',
        '-O', '3', '-C', '3', '-P', 'palette.txt', '-A', 'mean_j', '-g', 'annotation.gtf', '-F', 'pdf'
    ])
    if not (out.is_file() and out.stat().st_size > 0):
        failed.append(region)
assert not failed, f'no figure for {failed}'
```

MXE files also carry `upstreamES`/`downstreamEE`, and that span already covers both alternative exons. `examples/plot_sashimi.py` `batch_plot_rmats_events()` is the same recipe with the `--shrink` guard and a `RuntimeError` listing every failed event.

## rmats2sashimiplot

**Goal:** Plot directly from rMATS event coordinates without manual region calculation.

**Approach:** Filter the rMATS event file to the events to plot (it draws every row), pass BAM lists + a group file + event type, then check the output: rmats2sashimiplot **exits 0 when it fails** and leaves `Sashimi_plot/` empty. The contig is matched to the BAM header automatically (`chrX` in the event file, `X` in the BAM works).

```bash
# rmats2sashimiplot plots every row: keep only significant events (columns found by header name)
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)c[$i]=i; print; next}
    $c["FDR"]<0.05 && ($c["IncLevelDifference"]>0.1 || $c["IncLevelDifference"]<-0.1)' \
    rmats_output/SE.MATS.JC.txt > sig.SE.MATS.JC.txt

# group file: "label: first-last", 1-based over the --b1 replicates then the --b2 replicates
printf 'Control: 1-3\nTreatment: 4-6\n' > grouping.gf

rmats2sashimiplot \
    --b1 ctrl1.bam,ctrl2.bam,ctrl3.bam \
    --b2 trt1.bam,trt2.bam,trt3.bam \
    --event-type SE \
    -e sig.SE.MATS.JC.txt \
    --l1 Control \
    --l2 Treatment \
    -o sashimi_rmats \
    --exon_s 1 \
    --intron_s 5 \
    --group-info grouping.gf \
    --color '#1f77b4,#ff7f0e'

n_events=$(( $(wc -l < sig.SE.MATS.JC.txt) - 1 ))
n_pdf=$(find sashimi_rmats/Sashimi_plot -name '*.pdf' -size +0 2>/dev/null | wc -l)
[ "$n_pdf" -eq "$n_events" ] || { echo "rmats2sashimiplot wrote $n_pdf of $n_events figures" >&2; exit 1; }
```

`--event-type` (4.0.0; the old `-t SE` is rejected, rc 2) takes SE, A5SS, A3SS, MXE or RI. `--exon_s 1 --intron_s 5` draws introns at 1/5 of their real length. `--group-info` gives one plot per group (arc labels = group mean, plus the group's mean IncLevel); without it there is one plot per replicate and `--color` needs one colour per replicate, otherwise it prints `Error: Must provide sample label and color for each entry in bam_files!` and still exits 0.

## MAJIQ-VOILA Interactive Viewer

**Goal:** Browse LSV posterior PSI distributions interactively with splice-graph topology.

**Approach:** Run `voila view` on MAJIQ output; it starts a local web server (open the printed address in a browser; there is no `-o` output file).

MAJIQ/VOILA (bundled with MAJIQ, majiq.biociphers.org) is licence-gated (academic/commercial download) and was **not installed or run** in testing; the commands follow MAJIQ's public docs, so check `voila view --help` for your version.

```bash
# splicegraph file name and format depend on the MAJIQ version used for the build
voila view -p 5000 -j 8 build/splicegraph.<ext> psi_output/sample.psi.voila
voila view -p 5000 -j 8 build/splicegraph.<ext> deltapsi_output/group1_group2.deltapsi.voila
```

VOILA shows:
- Complete LSV graphs (single source / single target nodes)
- Per-junction posterior PSI violin plots
- ΔPSI distributions across all conditions
- Confidence by junction within an LSV

**The only tool that visualizes complex multi-junction LSVs intuitively.** For events that don't fit canonical SE/A5SS/A3SS, VOILA is the visualization of choice. It needs the MAJIQ build's splicegraph plus the `.voila` file; without a licence use ggsashimi on the region instead.

## leafviz Shiny App

**Goal:** Browse leafcutter clusters with intron-level effects and sashimi-like plots.

**Approach:** leafviz is a script directory inside the leafcutter repo (not an R package: `library(leafviz)` and `run_leafviz()` do not exist). Build annotation files from the GTF, prepare the results `.RData`, then launch the Shiny app from the `leafviz` directory. leafcutter is a GitHub R package (`devtools::install_github('davidaknowles/leafcutter/leafcutter')`), not Bioconductor (the as-shipped 0.2.9 fails to build against rstan >= 2.33 because of the old Stan array syntax; confirm `library(leafcutter)` loads); the `leafviz/` directory comes with its repo. Checked end to end on the planted 3v3 leafcutter results (leafcutter 0.2.9; the app serves HTTP 200).

```bash
# annotation_code = prefix of four files (_all_exons.txt.gz, _all_introns.bed.gz, _fiveprime.bed.gz, _threeprime.bed.gz);
# build them from the GTF version used in the differential analysis
perl leafcutter/leafviz/gtf2leafcutter.pl -o annot annotation.gtf

# groups.txt = the support file given to leafcutter_ds.R (sample <TAB> condition)
Rscript leafcutter/leafviz/prepare_results.R \
    -o leafviz.RData \
    -m groups.txt \
    leafcutter_perind_numers.counts.gz \
    ds_results_cluster_significance.txt \
    ds_results_effect_sizes.txt \
    annot

# runApp() uses the working directory: start from leafviz/, pass the .RData by absolute path
cd leafcutter/leafviz && Rscript run_leafviz.R /abs/path/leafviz.RData    # prints "Listening on http://127.0.0.1:<port>"
```

`download_human_annotation_codes.sh` in the same directory fetches prebuilt hg19 codes. Useful for cohort-level interactive filtering of clusters.

## Jutils for Tool-Agnostic Output

**Goal:** Visualize differential splicing output uniformly across rMATS, leafcutter, MntJULiP, and MAJIQ.

**Approach:** Convert tool output to Jutils' standard TSV, then plot. Run from the Jutils clone (`python3 jutils.py ...`). Run on rMATS output; the leafcutter/MntJULiP/MAJIQ converters follow `jutils.py convert-results --help` and were not run.

```bash
# writes rmats_JC_results.tsv and rmats_JCEC_results.tsv into --out-dir
python3 jutils.py convert-results --rmats-dir rmats_output/ --out-dir jutils_out/

# meta.tsv: sample<TAB>condition. Needs >= 2 events passing the cutoffs; writes clustermap*.pdf
python3 jutils.py heatmap --tsv-file jutils_out/rmats_JC_results.tsv --meta-file meta.tsv --q-value 0.05 --out-dir hm/ --pdf

# bam_list.tsv: sample<TAB>bam<TAB>condition
python3 jutils.py sashimi --tsv-file jutils_out/rmats_JC_results.tsv --meta-file meta.tsv \
    --gtf annotation.gtf --coordinate chr1:1000-2000 --bam-list bam_list.tsv --out-dir sh/ --pdf

# --tsv-file-list is a FILE with one "path<TAB>label" line per TSV, not a comma-separated list
printf 'jutils_out/rmats_JC_results.tsv\trMATS_JC\njutils_out/rmats_JCEC_results.tsv\trMATS_JCEC\n' > tsv_list.txt   # one line per TSV to compare
python3 jutils.py venn-diagram --tsv-file-list tsv_list.txt --out-dir vn/
```

(Yang 2021 *Bioinformatics*) Useful when comparing multiple tools' outputs across publications or doing meta-analysis. The sashimi labels are per-sample junction counts and matched pysam.

## pyGenomeTracks for Multi-Track Figures

**Goal:** Combine splicing with chromatin or coverage tracks for publication figures.

**Approach:** Build coverage bedGraphs and a junction BEDPE from the BAMs, define tracks in an INI file (genes, bedGraph/BigWig, BED, links), then run `pyGenomeTracks --tracks tracks.ini --region ... -o figure.pdf`. pyGenomeTracks 3.9 cannot draw a BAM (`InputError ... can not identify file type`).

```bash
# one merged BAM per group; -split is essential: without it introns are filled with coverage
samtools merge -f ctrl_merged.bam ctrl1.bam ctrl2.bam ctrl3.bam && samtools index ctrl_merged.bam
samtools merge -f trt_merged.bam trt1.bam trt2.bam trt3.bam && samtools index trt_merged.bam
bedtools genomecov -ibam ctrl_merged.bam -split -bga > ctrl.bedgraph
bedtools genomecov -ibam trt_merged.bam -split -bga > trt.bedgraph
```

```ini
[gene_models]
file = annotation.gtf
height = 3
title = GENCODE v45
fontsize = 10
file_type = gtf

[ctrl_coverage]
file = ctrl.bedgraph
title = Control
color = #1f77b4
height = 3
min_value = 0
max_value = 200
file_type = bedgraph

[trt_coverage]
file = trt.bedgraph
title = Treatment
color = #ff7f0e
height = 3
min_value = 0
max_value = 200
file_type = bedgraph

[junctions]
file = junctions.bedpe
title = Junctions
height = 2
file_type = links
links_type = arcs
```

Tracks are scaled independently: set the same `min_value`/`max_value` on both coverage tracks (pick `max_value` from the data) or the two groups are not comparable. A BigWig made from the same `-split` bedGraph works too (`file_type = bigwig`).

The `junctions.bedpe` file must be in **BEDPE format** (6 columns: chr1 start1 end1 chr2 start2 end2 [+ optional score]). Convert from regtools .bed12 junctions (the score is the read count summed over the merged BAM; `-s XS` needs XS-tagged BAMs, otherwise the strand is `?`):

```bash
samtools merge -f all_merged.bam ctrl_merged.bam trt_merged.bam && samtools index all_merged.bam   # regtools needs an indexed BAM
regtools junctions extract -s XS -o regtools_junctions.bed all_merged.bam
# regtools BED12 column 11 is blockSizes (anchor_left, anchor_right);
# column 12 is blockStarts (0, intron_length + anchor_left).
# Intron start = chromStart + anchor_left = $2 + a[1]
# Intron end   = chromStart + blockStarts[2] = $2 + b[2]
awk 'BEGIN{OFS="\t"} {split($11,a,","); split($12,b,","); s=$2+a[1]; e=$2+b[2]; print $1, s, s+1, $1, e-1, e, $5}' \
    regtools_junctions.bed > junctions.bedpe
```

```bash
pyGenomeTracks --tracks tracks.ini --region chr17:43094000-43125000 -o figure.pdf
```

## Reading Sashimi Plots (Interpretation Guide)

| Visual element | What it represents |
|----------------|--------------------|
| Filled coverage track | Read coverage at each genomic position (per sample, overlaid; the group mean with `-A mean`) |
| Arc / curve between exons | Junction-spanning reads; arc connects donor to acceptor |
| Number on arc | Count of junction-spanning reads (raw per sample; with `-A` the rounded group mean) |
| Arc thickness | Often proportional to read count (tool-dependent) |
| Gene model below | Exons (boxes) and introns (lines) from GTF |
| Multiple parallel tracks | Per-sample (default) or per-group (with `-O`) |

**Junction count interpretation:** the number on an arc is the absolute count of reads whose CIGAR string contained an `N` operation matching that intron coordinate. Higher = more usage. Compare counts on inclusion vs skipping arcs to estimate PSI visually.

## Per-Tool Failure Modes

### ggsashimi: Off-Strand Junction Artifacts

**Trigger:** Stranded RNA-seq library plotted without strand specification.

**Mechanism:** ggsashimi reads BAM strand from CIGAR + flag; without strand info, antisense junctions appear as artifacts.

**Symptom:** Implausible junctions in regions with overlapping antisense genes; "noise" arcs at unexpected locations.

**Fix:** Set library strandedness with `-s MATE2_SENSE` (dUTP/TruSeq reverse-stranded PE; `-s MATE1_SENSE` for forward PE; verify orientation with RSeQC `infer_experiment.py`). `MATE1_SENSE`/`MATE2_SENSE` are paired-end only (`TypeError: ... 'NoneType' and 'bool'` on single-end); single-end uses `-s SENSE`/`ANTISENSE`. With `-s` the output is two files, `<prefix>_+.<fmt>` and `<prefix>_-.<fmt>`. Alternatively, pre-filter the BAM by strand with `samtools view -f 16` / `-F 16`.

### ggsashimi: Silent Failures (exit 0)

**Trigger:** BAM path typo in the TSV, region without reads, R package missing or ggplot2 error (e.g. `Unknown colour name`).

**Symptom:** rc 0 with a dropped sample, an empty figure, or no figure at all.

**Fix:** Check that every BAM exists, that the region has reads (`samtools view -c sample.bam chr1:100-200`), and that the figure file exists and is non-empty (recipes above).

### leafviz: Annotation Codes Mismatch

**Trigger:** Using leafviz with annotation_codes from different GENCODE version than leafcutter clusters.

**Mechanism:** annotation_codes encodes intron-to-event-class mapping per GTF version.

**Symptom:** Many clusters show as "unannotated" despite being in canonical GTF.

**Fix:** Generate annotation_codes with `gtf2leafcutter.pl` from the same GTF used in differential analysis.

## Best Practices

| Tip | Rationale |
|-----|-----------|
| Use `--shrink` for genes with large introns | Keeps exons visible (TTN, brain genes with multi-kb introns) |
| Limit to 3-4 groups per figure | More becomes hard to read |
| Include 200-500 nt flanking exons | Show full splicing context |
| For MXE events, plot both alternative exons | Otherwise only half of the event is visible |
| Check accessibility colors | Use ColorBrewer-safe palettes for color-blind readers |
| Always include a legend | Sashimi figures without legends are uninformative for non-experts |

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| R error printed (`there is no package called ...`, `Unknown colour name`), **rc 0, no figure** | ggsashimi needs `R` on PATH with ggplot2/data.table/gridExtra/gtable (the samtools binary is not used; pysam is); or a bad palette colour | Install the R packages; assert the figure exists |
| `ggsashimi: ValueError: invalid contig` | Contig name differs between region and BAM (`chr1` vs `1`) | Map the name against the BAM header (recipes above) |
| `ERROR: No available bam files.` | Every path in the TSV is wrong (a single wrong path is dropped silently) | Check paths; relative paths resolve against the TSV's directory |
| `ERROR: Cannot apply aggregate function if overlay is not selected.` | `-A` without `-O` | Add `-O 3` |
| `RuntimeError: generator raised StopIteration` | `--shrink` with no junction passing `-M` | Lower `-M` or drop `--shrink` |
| `rmats2sashimiplot: unrecognized arguments: -t SE` | Old flag | `--event-type SE` |
| `Error: Must provide sample label and color for each entry in bam_files!`, rc 0, empty `Sashimi_plot/` | Colours per replicate given for 2 groups | Add `--group-info grouping.gf` or one colour per replicate |
| `pyGenomeTracks: InputError ... can not identify file type` | BAM track (unsupported) or missing `file_type` | Use bedGraph/BigWig; set `file_type` |
| `jutils.py venn-diagram: FileNotFoundError` on `a.tsv,b.tsv` | `--tsv-file-list` is a file of paths | Write the list file |
| `App dir must contain either app.R or server.R` | `run_leafviz.R` not started from the `leafviz/` directory | `cd leafcutter/leafviz` first |
| `prepare_results.R`: `<file> does not exist` | annotation_code prefix or a path is wrong | Re-run with correct paths |

## Troubleshooting

| Issue | Cause | Solution |
|-------|-------|----------|
| No junctions shown | `-M` above every junction's count (default is 1) | Lower `-M` |
| Plot too crowded | Many samples without aggregation | Use `-O 3` to overlay groups |
| Annotation missing or wrong gene | GTF lacks gene_name attribute or wrong build | Verify GTF version vs BAM reference; pre-filter the GTF to the relevant features |
| Memory issues on large regions | >100 kb regions with many samples | Plot smaller windows or pre-extract reads with samtools view |
| Y-axis dominated by one peak | Outlier sample | Filter the outlier out of the TSV |
| Gene model shifted against coverage, tick labels clipped | ggplot2 >= 3.5 | Use ggplot2 3.4.4 |

## Related Skills

- differential-splicing - Identify events to plot; sashimi plots are validation
- splicing-quantification - Context for PSI values; sashimi provides visual confirmation
- data-visualization/genome-tracks - Multi-track figure design (pyGenomeTracks, Gviz)
- data-visualization/ggplot2-fundamentals - ggsashimi customization (extends ggplot2)
- data-visualization/color-palettes - Accessible color choices
- data-visualization/volcano-and-ma-plots - Volcano complement to sashimi
- data-visualization/heatmaps-clustering - Heatmap complement to sashimi

## References

- Katz et al 2010 *Nat Methods* - MISO sashimi plot original
- Garrido-Martin et al 2018 *PLoS Comput Biol* - ggsashimi
- Yang et al 2021 *Bioinformatics* - Jutils
- Vaquero-Garcia et al 2016 *eLife* - MAJIQ / VOILA
- Li et al 2018 *Nat Genet* - leafcutter / leafviz
- Ramirez et al 2018 *Nat Commun* - pyGenomeTracks
