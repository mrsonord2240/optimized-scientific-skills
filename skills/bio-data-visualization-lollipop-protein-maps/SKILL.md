---
name: bio-data-visualization-lollipop-protein-maps
description: Plot per-gene mutation distributions on protein-domain maps with maftools, trackViewer, and g3viz. Use for recurrence-aware mutation maps, custom UniProt domains, or an interactive HTML supplement.
tool_type: mixed
primary_tool: maftools
license: MIT
author: GPTomics
---

# Lollipop / Needle Protein Maps

Use a lollipop plot to show where recurrent protein-altering calls fall along one protein. In maftools, stem height is the number of MAF rows for one protein-change string; it is not population frequency, point size, or a deduplicated sample count. State the chosen recurrence unit in the caption and validate a novel hotspot in an independent cohort before interpreting it.

## Choose a route

| Need | Use | Notes |
| --- | --- | --- |
| Fast MAF plot or paired cohorts | `maftools::lollipopPlot` / `lollipopPlot2` | Base graphics; returns a data table, not a ggplot object. |
| Exact, current domain coordinates or custom labels | `trackViewer::lolliplot` | Supply a protein-coordinate `GRanges`; use UniProt features. |
| Interactive supplement | `g3viz::g3Lollipop` | Write the widget with `htmlwidgets::saveWidget`. g3viz 1.2.0 is archived on CRAN. |

## Input and isoform checks

1. Keep the MAF's protein-change column and check it before plotting. Do not force `AACol='HGVSp_Short'` unless that column exists; maftools can auto-detect supported `HGVSp_Short`, `Protein_Change`, or `AAChange` columns.
2. Record the transcript/protein used by the variant caller and the plot. maftools' bundled domain table selects the longer listed RefSeq transcript when no ID is supplied; it is not a live Pfam query or a promise of the canonical UniProt isoform. Use a valid RefSeq mRNA ID when a particular transcript is required, e.g. TP53 `refSeqID='NM_000546'`. `proteinID` takes an NCBI protein ID such as `NP_000537`, not a UniProt accession such as `P04637`.
3. TP53-201's Ensembl transcript is `ENST00000269305`. Do not use `ENST00000288602`: that is BRAF-201.
4. If calls extend beyond a requested shorter isoform, maftools can clip them without a warning. Compare the largest parsed position with the plotted protein length first.

```r
library(maftools)
maf <- read.maf("cohort.maf", verbose = FALSE)
change_col <- intersect(c("HGVSp_Short", "Protein_Change", "AAChange"), names(maf@data))
if (!length(change_col)) stop("No supported protein-change column in MAF")
change_col <- change_col[[1L]]

pdf("TP53_lollipop.pdf", width = 8, height = 4)
result <- lollipopPlot(
  maf = maf, gene = "TP53", AACol = change_col, refSeqID = "NM_000546",
  labelPos = c(175, 248, 273), labPosSize = 1,
  showMutationRate = TRUE, domainLabelSize = 0.8,
  printCount = TRUE, # prints the per-change table to the console; it does not label the figure
  colors = c(Missense_Mutation = "#D55E00", Nonsense_Mutation = "#000000",
             Frame_Shift_Del = "#0072B2", Frame_Shift_Ins = "#56B4E9",
             Splice_Site = "#CC79A7", In_Frame_Del = "#009E73",
             In_Frame_Ins = "#F0E442", Translation_Start_Site = "#999999",
             Nonstop_Mutation = "#E69F00")
)
dev.off()
```

`printCount=TRUE` prints a table for reporting; use `trackViewer` (below) or base-graphics `text()` after inspecting the drawn coordinates when a figure must carry count labels. maftools draws one stem per distinct protein change: KRAS G12D and G12V remain separate, rather than becoming a single 55-call G12 stem. Deduplicate a duplicate call on `Tumor_Sample_Barcode + HGVSp` before plotting if the intended height is unique samples.

```r
maf_a <- subsetMaf(maf, clinQuery = 'Subtype == "Luminal"')
maf_b <- subsetMaf(maf, clinQuery = 'Subtype == "Basal"')
pdf("TP53_lollipop_subtype.pdf", width = 10, height = 5)
lollipopPlot2(m1 = maf_a, m2 = maf_b, gene = "TP53", m1_name = "Luminal", m2_name = "Basal",
              AACol1 = change_col, AACol2 = change_col)
dev.off()
```

For crowded plots, reduce `domainLabelSize`, use a few `labelPos` values, and check the rendered output. Domain colours in maftools' bundled table are not a functional ontology and may vary between calls; use trackViewer when stable functional colours matter.

## Parse protein positions safely

HGVSp strings are heterogeneous. The helper below accepts one-letter (`p.R175H`) and three-letter (`p.Arg175His`) forms plus stop-loss extensions (`p.*394Wext*?`), returns `NA` for silent/unknown/empty strings, and warns with the exact dropped rows. Do not silently feed `NA` positions to `IRanges`.

```r
protein_position <- function(change) {
  x <- trimws(as.character(change))
  take <- function(pattern, x) {
    hits <- regmatches(x, regexec(pattern, x))
    vapply(hits, function(h) if (length(h) >= 2L) as.integer(h[[2L]]) else NA_integer_, integer(1))
  }
  out <- take("^p\\.[A-Z*](\\d+)", x)
  out[grepl("^p\\.[A-Z][0-9]+\\?$", x)] <- NA_integer_ # p.M1? is uncertain, unlike p.*394Wext*?
  missing <- is.na(out)
  out[missing] <- take("^p\\.[A-Z][a-z]{2}(\\d+)", x[missing])
  out
}
```

## trackViewer with inclusive protein coordinates

Protein feature coordinates are inclusive. `IRanges(start, end)` therefore preserves both endpoints; do not calculate a width by hand. The supplied TP53 P04637 map combines current UniProt feature ranges for Transactivation 1--44, oligomerization 325--356, and the basic region 368--387 with the conventional core DNA-binding range 102--292. Current UniProt P04637 JSON does not label 102--292 as a single DNA-binding feature, so describe it as a supplied conventional range and cite the chosen literature source if that boundary is used in a publication. Confirm the accession and feature release as well.

```r
library(data.table)
library(trackViewer)
library(GenomicRanges)

class_col <- c(Missense_Mutation = "#D55E00", Nonsense_Mutation = "#000000",
               Frame_Shift_Del = "#0072B2", Frame_Shift_Ins = "#56B4E9",
               Splice_Site = "#CC79A7", In_Frame_Del = "#009E73",
               In_Frame_Ins = "#F0E442", Translation_Start_Site = "#999999",
               Nonstop_Mutation = "#E69F00")
calls <- as.data.table(maf@data)[Hugo_Symbol == "TP53"]
calls[, aa_pos := protein_position(get(change_col))]
bad <- calls[is.na(aa_pos)]
if (nrow(bad)) warning("Dropping ", nrow(bad), " TP53 call(s) without a parseable protein position: ", paste(unique(bad[[change_col]]), collapse = ", "))
calls <- calls[!is.na(aa_pos) & aa_pos >= 1L & aa_pos <= 393L]
unknown <- setdiff(unique(as.character(calls$Variant_Classification)), names(class_col))
if (length(unknown)) warning("Using grey for unmapped variant class(es): ", paste(unknown, collapse = ", "))
calls[, plot_class := fifelse(as.character(Variant_Classification) %in% names(class_col), as.character(Variant_Classification), "Other")]
class_col <- c(class_col, Other = "#666666")
summary <- calls[, .(count = .N, class = plot_class[1L], residue = unique(sub("^p\\.", "", get(change_col)))[1L]), by = aa_pos]
snps <- GRanges("TP53", IRanges(summary$aa_pos, width = 1L), color = unname(class_col[as.character(summary$class)]), score = summary$count)
names(snps) <- ifelse(summary$count >= 2L, summary$residue, "")
# The 102--292 DNA-binding core is a conventional supplied range; the other
# three ranges above are current UniProt P04637 features.
features <- GRanges("TP53", IRanges(c(1, 102, 325, 368), c(44, 292, 356, 387), names = c("Transactivation", "DNA binding core (conventional)", "Oligomerization", "Basic")), fill = c("#56B4E9", "#0072B2", "#009E73", "#CC79A7"), height = 0.04)
pdf("TP53_trackviewer.pdf", width = 10, height = 4)
lolliplot(snps, features, ylab = "Mutation-row count", xaxis = TRUE, yaxis = TRUE, legend = list(labels = names(class_col), col = unname(class_col)))
dev.off()
```

The summary intentionally aggregates by residue only for this custom plot. If several classes occur at one residue, choose and document a rule (for example, the first observed class, a multi-class annotation, or separate tracks).

## Interactive HTML with g3viz

`g3viz` 1.2.0 was archived from CRAN; an environment that already supplies it can use the following route. `install.packages("g3viz")` may fail; install a reviewed CRAN Archive source only under the environment's package policy. `output.filename` names the widget's download button, not an HTML output path.

```r
library(g3viz)
library(htmlwidgets)
mut <- readMAF("cohort.maf", protein.change.col = change_col)
widget <- g3Lollipop(mut, gene.symbol = "TP53", protein.change.col = change_col, plot.options = g3Lollipop.theme(theme.name = "nature"), output.filename = "TP53_lollipop")
saveWidget(widget, "TP53_lollipop.html", selfcontained = TRUE)
```

## Interpretation and failure checks

| Check | Why it matters | Action |
| --- | --- | --- |
| Rows versus samples | Duplicate calls inflate a stem. | Deduplicate sample + change for sample recurrence; report the denominator. |
| One change versus one residue | maftools keeps G12D and G12V as separate stems. | Aggregate deliberately only when residue recurrence is the question. |
| `p.M1?`, `p.=`, empty HGVSp | These have no usable amino-acid position. | Warn and report the dropped count. |
| 3-letter or stop-extension HGVSp | Simplistic one-letter regex returns `NA`. | Parse explicitly and test the result. |
| Isoform length | Positions past a shorter selected protein can disappear. | Compare max position with length and state the ID. |
| Novel hotspot | Coverage and cohort composition can mimic recurrence. | Replicate in an independent cohort; use a formal hotspot method for inference. |

## References

- Mayakonda A, et al. 2018. *Genome Research* 28:1747-1756 (maftools).
- Ou J, Zhu LJ. 2019. *Nature Methods* 16:453-454 (trackViewer).
- Chang MT, et al. 2016. *Nature Biotechnology* 34:155-163 (recurrent mutations).
- Lawrence MS, et al. 2014. *Nature* 505:495-501 (hotspot inference).

## Related Skills

- data-visualization/oncoprint-mutation-matrices
- variant-calling/variant-annotation
- data-visualization/color-palettes
