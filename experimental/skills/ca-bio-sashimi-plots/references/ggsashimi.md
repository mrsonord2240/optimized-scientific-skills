# ggsashimi options and checks


**Goal:** Generate publication-quality sashimi plot for a region with samples grouped by condition and per-sample tracks aggregated.

**Approach:** Define samples + groups in a TSV (no header), a palette file, then call ggsashimi with coordinates, GTF, and visual flags. ggsashimi exits 0 on several failures (see Silent Failures): `plot_sashimi()` checks the inputs before and the figure after; on the bare command line, do it yourself.

```python
# scripts/plot_sashimi.py: create_grouping_file, write_palette, plot_sashimi (run from examples/, or put it on sys.path)
from plot_sashimi import create_grouping_file, write_palette, plot_sashimi

bams = ['ctrl1.bam', 'ctrl2.bam', 'ctrl3.bam', 'trt1.bam', 'trt2.bam', 'trt3.bam']
create_grouping_file(bams, ['Control'] * 3 + ['Treatment'] * 3, 'sashimi_groups.tsv')  # col1 sample id, col2 BAM, col3 group (-O overlay, -C colour)
write_palette(['#1f77b4', '#ff7f0e'], 'palette.txt')   # one colour per group, in order of first appearance
plot_sashimi('sashimi_groups.tsv', 'chr17:43094000-43125000', 'BRCA1_sashimi', 'gencode_v45.gtf',
             options={'palette': 'palette.txt', 'width': 10})   # contig spelled as in the BAM header (or mapped for you)
```

The equivalent bare command line, which also reaches `--ann-height` and `--base-size` (not exposed by `plot_sashimi()`); the contig must be spelled as in the BAM header:

```bash
ggsashimi.py -b sashimi_groups.tsv -c chr17:43094000-43125000 -o BRCA1_sashimi -g gencode_v45.gtf \
    --alpha 0.25 --height 3 --width 10 --ann-height 4 --base-size 14 --shrink --fix-y-scale \
    -O 3 -C 3 -P palette.txt -A mean_j -F pdf
test -s BRCA1_sashimi.pdf || echo 'no figure written (R error above?)' >&2
```

`plot_sashimi()` adds to the command line: a missing-BAM check, contig mapping (chrX vs X), an empty-region check, the `--shrink` guard below and a figure-exists check. Checked: ggsashimi's junction labels equal an independent pysam count (planted 3v3 set, real ENCODE 12-BAM locus, real chrX BAMs).

Key ggsashimi flags (Garrido-Martin 2018 *PLoS Comput Biol*):
- `-O 3`: column 3 of the TSV is the overlay level; samples of a group are drawn in one track. Required for `-A`
- `-C 3 -P palette.txt`: colour by column 3 using the palette file (R colour names or hex, one per line). Without `-C` everything is grey; with `-C` and no `-P` the colours are R defaults (red/green), not blue/orange
- `-A mean_j`: the arc label is the rounded (half to even: 6.5 shows 6) plain mean of the raw junction counts of the group's samples that have the junction; every sample's coverage stays overlaid. `mean` and `median` also aggregate the coverage; `median_j` is the median analogue. There is no depth normalization
- `-M`: minimum reads for a junction to be drawn, **default 1**, inclusive, applied **per sample before `-A`**: samples below `-M` drop out of the mean, biasing labels upward (planted 10.33 shows 11 at `-M 10`; 10 at `-M 1`). Keep `-M 1` for exact labels; raise it (5-10, 20+ for crowded plots) only to declutter. If `-M` is above every junction's count, ggsashimi still exits 0 and writes a coverage-only figure with no arcs (`scripts/plot_sashimi.py` warns)
- `--shrink`: rescale long introns for compact display (keeps exons visible in genes with multi-kb introns). Crashes with `RuntimeError: generator raised StopIteration` when no junction passes `-M` (with `-s`, when one strand has none): lower `-M` or drop `--shrink`
- `--fix-y-scale`: identical y-axis across groups (essential for visual comparison)
- `--alpha 0.25`: transparency of per-sample coverage in overlay mode
- `--height`/`--width`/`--ann-height`/`--base-size`: sizes in inches / font size (e.g. `--width 12 --height 4 --base-size 16`)
- `-F pdf|svg|png|jpeg|tiff` (`-R` for raster PPI); pick explicitly: PDF for publication, PNG for slides, SVG for editing
- `-g`: GTF with exons; ggsashimi has no feature filter, so pre-filter the GTF (`awk '$0 ~ /protein_coding/'`) to restrict it
- Colour convention: control = blue (`#1f77b4`), treatment = orange (`#ff7f0e`); document it, use ColorBrewer for >2 groups

## Batch Plotting from rMATS Hits

**Goal:** Auto-generate sashimi plots for all significant rMATS differential events.

**Approach:** `batch_plot_rmats_events()` parses SE.MATS.JC.txt, keeps the top events by significance and effect size, expands each to `upstreamES`-500 .. `downstreamEE`+500, maps the contig name onto the BAM header, clamps the start to 1, runs `plot_sashimi()` per event (file names sanitised, event ID included) and raises a `RuntimeError` listing every event that produced no figure. rMATS writes `chrX`; an Ensembl-style BAM calls it `X` and ggsashimi dies with `ValueError: invalid contig`. Events near a contig start give a start < 1 (handled).

```python
from plot_sashimi import batch_plot_rmats_events    # same import as above; reuses sashimi_groups.tsv and palette.txt
batch_plot_rmats_events('rmats_output/SE.MATS.JC.txt', 'sashimi_groups.tsv', 'annotation.gtf', 'sashimi_plots/',
                        n_top=25, fdr_cutoff=0.05, dpsi_cutoff=0.1, flank=500, palette='palette.txt')
```

MXE files also carry `upstreamES`/`downstreamEE`, and that span already covers both alternative exons. `--shrink` is dropped with a warning for an event with no junction at `-M` (see the flag list).
