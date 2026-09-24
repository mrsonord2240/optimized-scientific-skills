# Sequence Logos - Usage Guide

## Overview

Sequence logos visualize per-position base/aa composition for an aligned motif set. Position height encodes information content in bits (Schneider-Stephens 1990); letter height is proportional to frequency. The bits encoding is the canonical bioinformatics standard - conserved positions are visually tall, variable positions short. ggseqlogo (R) is the modern default; Logomaker (Python) is the most flexible for non-standard alphabets; WebLogo is the CLI / web original.

## Prerequisites

```r
install.packages(c('ggseqlogo', 'patchwork'))
```

```bash
pip install logomaker
# CLI option:
pip install weblogo
# On Windows, install Ghostscript for WebLogo PDF/PNG/SVG output.
```

## Quick Start

Tell your AI agent what you want to do:
- "Plot a sequence logo of these aligned TF binding sites with bits encoding"
- "Apply human-genome background composition with a relative-entropy logo"
- "Stack two motifs vertically to compare TF-A vs TF-B"
- "Custom color scheme for a protein motif by amino-acid class"
- "Show enrichment AND depletion using log-odds encoding"

## Example Prompts

### TF binding-site logo

> "Build a relative-entropy logo of CTCF binding sites against human genome composition (A=0.29, C=0.21, G=0.21, T=0.29). Use a ggseqlogo custom matrix and label the y-axis as bits relative to background."

### Splice site composition

> "Logo of 5' splice sites from 200 aligned exon-intron junctions. Bits encoding, uniform background. ggseqlogo."

### Differential motif

> "Use Logomaker weight (log-odds) encoding to show enrichment AND depletion in TF-A motif relative to TF-B background. Flip depleted letters below the axis."

### Protein kinase substrate

> "Logo of phosphorylation-site neighborhoods (15 aa around phospho-S/T). Color by amino-acid class: phospho-acceptor red, basic blue, acidic purple."

### Multi-motif comparison

> "Stack three TF logos (CTCF, REST, GATA1) vertically to compare conservation patterns. Same scale for all."

## What the Agent Will Do

1. Load aligned sequences (FASTA) or PWM (counts/probability/PSSM).
2. Verify alignment: same length for all entries; sequence type (DNA, RNA, protein).
3. Choose uniform-background information or genome-background relative entropy.
4. Compute the stated per-position heights before plotting.
5. Render with ggseqlogo / Logomaker / WebLogo per project preference.
6. Annotate N (number of input sequences) for transparency.
7. Apply CVD-safe color scheme; for proteins use functional-class coloring.

## Tips

- **Default to bits**, not probability. Bits show the conservation gradient; probability shows raw frequency with every position equal-height.

- **Background matters.** ggseqlogo does not have a background-frequency argument. For its non-uniform-background route, compute raw, unsmoothed maximum-likelihood relative-entropy heights and use method = 'custom' (see examples/relative_entropy_logo.R). Logomaker accepts background; WebLogo accepts --composition. Choose and report any pseudocount or biological prior explicitly.

- **N ≥ 20 for credibility.** Below this, even random sequences look conserved. Record the tool's pseudocount/prior or correction, and treat N < 5 as exploratory rather than authoritative.

- **PWM orientation matters.** ggseqlogo expects letters as rows and positions as columns; Logomaker expects positions as rows and letters as columns. Verify dimensions before plotting.

- **Signed weight logos are not EDLogo.** Logomaker can display signed log-odds below the axis, but use an actual EDLogo implementation if its specific enrichment/depletion methodology is required.

- **Stacking logos requires same alphabet and same scale.** DNA (max 2 bits) vs protein (max 4.3 bits) cannot be visually compared without normalization.

- **Custom alphabets.** ggseqlogo auto-detects ordinary RNA U; set seq_type = 'rna' for clarity and use --sequence-type rna in WebLogo. Define explicit schemes for extended alphabets.

- **For PSSM input** (signed log-odds), pass a correctly constructed signed matrix to Logomaker; do not pretend an already weighted PSSM is counts.

- **WebLogo is best for batch CLI use** with `--composition equiprobable` (uniform) or `--composition <species>` for built-in genome compositions.

## Related Skills

- chip-seq/motif-analysis - Motif discovery that produces the PWM
- atac-seq/footprinting - Footprinting motifs to visualize
- clip-seq/clip-motif-analysis - CLIP-derived motifs
- alignment/multiple-alignment - Aligned sequences as input
- data-visualization/color-palettes - Custom alphabet schemes
