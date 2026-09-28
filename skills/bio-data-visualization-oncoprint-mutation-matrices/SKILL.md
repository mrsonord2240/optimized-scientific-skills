---
name: bio-data-visualization-oncoprint-mutation-matrices
description: Build OncoPrint and co-mutation matrix plots from somatic-variant cohorts using ComplexHeatmap, maftools, and comut.py with alteration-type stacking, sample ordering by mutational burden, mutual-exclusivity overlays, and clinical annotation tracks. Use when visualizing per-sample mutation patterns across recurrent driver genes, comparing alteration classes, or identifying mutually-exclusive / co-occurring driver pairs.
tool_type: mixed
primary_tool: ComplexHeatmap
---

## Version Compatibility

The shipped paths were checked with ComplexHeatmap 2.22.0, maftools 2.22.0,
circlize 0.4.18, comut 0.0.3, pandas 2.3.3, and matplotlib 3.10. comut 0.0.3's
continuous track is incompatible with pandas 3; run `scripts/comut_plot.py` in a
pandas 2.x environment.

Before adapting a path, check `packageVersion()` in R or `pip show` in Python.
Install missing R packages with
`BiocManager::install(c("ComplexHeatmap", "maftools"))` plus
`install.packages("circlize")`; install the Python path with
`pip install "pandas<3" comut matplotlib`.

# OncoPrint and Mutation Matrix Plots

An OncoPrint is a gene-by-sample matrix in which a cell can contain several
alteration classes. Preserve those classes: a copy gain plus a missense event is
not one categorical state.

## Choose the workflow and ordering

| Question | Tool and ordering | Display |
|---|---|---|
| Which genes are most altered? | ComplexHeatmap with `row_order = order(-rowSums(mat != ''))` | Cohort-wide gene percentages and sample burden |
| Canonical mutation staircase | ComplexHeatmap default column memo sort | Binary alteration pattern across top genes |
| Per-patient burden | Explicit `column_order = order(-clinical$tmb)` | TMB track and sample labels when readable |
| Subtype-driver enrichment | `column_split = clinical$Subtype` after exact ID alignment | Split columns; percentages and right bar remain cohort-wide |
| Rapid TCGA MAF view | maftools `oncoplot()` | Only samples represented in the MAF |
| Python workflow | `scripts/comut_plot.py` | Burden-sorted, cohort-complete comut plot |
| Co-occurrence or mutual exclusion | maftools `somaticInteractions()` | Pairwise Fisher results plus the plotted signed significance |

ComplexHeatmap's default row order counts alteration types, not unique mutated
samples. A multi-class cell can therefore move a less-prevalent gene above a
more-prevalent one. Pass the explicit `row_order` above when “most altered” means
number of samples.

## Build a cohort-complete matrix first

`scripts/maf_to_oncoprint.R` is the canonical preparation seam. It maps coding
MAF consequences, collapses duplicate same-class calls, accepts optional
normalized Amp/HomDel/Fusion calls, includes zero-mutation samples from an
explicit cohort list, rejects an all-empty selection, and aligns clinical rows
by exact sample ID.

```r
source("scripts/maf_to_oncoprint.R")
maf_df <- read.delim("cohort.maf", comment.char = "#", check.names = FALSE)
clinical <- read.delim("clinical.tsv", check.names = FALSE)
mat <- maf_to_oncoprint(maf_df, clinical$Tumor_Sample_Barcode, top = 20)
clinical <- align_oncoprint_clinical(clinical, colnames(mat))
stopifnot(identical(clinical$Tumor_Sample_Barcode, colnames(mat)))
```

The default map is explicit:

- `Missense_Mutation`, `In_Frame_Ins`, `In_Frame_Del` -> `Missense`
- `Nonsense_Mutation`, `Frame_Shift_Ins`, `Frame_Shift_Del`,
  `Nonstop_Mutation`, `Translation_Start_Site` -> `Truncating`
- `Splice_Site` -> `Splice`
- `Silent`, intronic, RNA, IGR, flank, and other unmapped classifications are
  excluded and attached to the matrix as `ignored_variant_classifications`.
- Supply CNV/fusion rows through `additional_calls` with columns
  `Hugo_Symbol`, `Tumor_Sample_Barcode`, and `alteration_class`.

Do not rely on row position for annotations. ComplexHeatmap accepts positional
vectors without warning, so a shuffled clinical table can silently label the
wrong samples.

## ComplexHeatmap rendering

Run the complete cohort-aware implementation as:

```bash
Rscript examples/oncoprint_phd.R cohort.maf clinical.tsv oncoprint.pdf
```

The example defines one renderer per class. Amp and HomDel fill the cell;
Missense, Truncating, and Splice use distinguishable partial-height rectangles;
Fusion is a triangle. It retains all cohort columns, log-transforms TMB, aligns
annotations, and orders rows by mutated-sample frequency.

ComplexHeatmap percentages are calculated from every column in the input
matrix. `remove_empty_columns = TRUE` hides empty displayed columns but does not
recompute those percentages. Use `FALSE` to keep the cohort visible; construct
`mat` from the explicit cohort list to make its denominator the intended cohort.

## maftools quick path

`clinicalData` must contain `Tumor_Sample_Barcode`. Give every discrete level a
color; partial maps silently turn unlisted groups grey. Numeric annotations need
a sequential palette.

```r
library(maftools)
maf <- read.maf("cohort.maf", clinicalData = clinical)
oncoplot(maf, top = 20,
         clinicalFeatures = c("Subtype", "Stage"),
         annotationColor = list(
           Subtype = c(Luminal="#0072B2", Basal="#D55E00", HER2="#009E73"),
           Stage = c(I="#FFFFCC", II="#FED976", III="#FD8D3C", IV="#BD0026")),
         sortByAnnotation = TRUE, removeNonMutated = FALSE)
```

maftools only knows samples present in the MAF. `removeNonMutated = FALSE` does
not add cohort members with no MAF row, even if they occur in `clinicalData`.
Use the explicit matrix plus ComplexHeatmap when the full cohort is the required
denominator. maftools also uses black for `Multi_Hit`; do not equate that black
tile with this Skill's black `Truncating` tile when comparing tools.

## comut.py path

Prepare tab-separated mutation, clinical, and TMB files with columns
`sample`, `category`, `value`, plus a cohort file with `sample`. The script sets
the full cohort before adding datasets, sorts samples by TMB (or mutation burden
without TMB), reverses the y-axis category order so the most frequent gene is on
top, derives the TMB range from the data with a zero lower bound, and adds a
unified alteration/clinical legend. Cohort and track IDs are normalized to
stripped strings; empty or NA IDs fail. TMB values must be finite and
non-negative.

```bash
python scripts/comut_plot.py --mutations mutations.tsv --cohort cohort.tsv --clinical clinical.tsv --tmb tmb.tsv --output comut.pdf
```

Sample IDs are shown for cohorts of at most 50 samples and hidden automatically
above that threshold. Change the cutoff with `--sample-label-threshold N`; use
`0` to hide all labels or a larger value to force labels for a known-readable
layout.

Unknown cohort IDs and alteration classes fail loudly. Use pandas 2.x; do not
work around the pandas 3 error by dropping the continuous track silently.

## Mutual exclusivity and co-occurrence

```r
si <- somaticInteractions(maf, top = 20,
                          pvalue = c(0.05, 0.01), fontSize = 0.7)
si[, c("gene1", "gene2", "pValue", "oddsRatio", "Event", "pAdj")]
```

The return value is a `data.table` with pair labels, p-values, odds ratios,
2-by-2 cell counts (`00`, `01`, `11`, `10`), adjusted p-values, and `Event`.
The signed `-log10(p)` representation belongs to the plot; it is not the return
object.

Fisher tests ignore sample-specific mutation-rate background. Prefer DISCOVER
for pan-cancer analyses with large burden heterogeneity. At N < 50, treat the
plot as descriptive, report per-gene exact-binomial intervals, and do not claim
pairwise mutual exclusion from sparse cells. A 0.5 continuity correction can
even reverse the apparent direction in very sparse tables; report raw 2-by-2
counts and avoid classifying direction from a corrected odds ratio alone.

## Quantitative guidance

| Item | Guidance |
|---|---|
| Displayed genes | Usually 10-25 in one panel |
| Sample labels | Hide when the cohort is too dense to read |
| Hypermutator track | `log10(tmb + 1)` or a clearly disclosed cap |
| Pair testing | Prefer N >= 100; use effect-size focus from 50-99 |
| Very small cohorts | Exact frequency intervals; no discovery claim from mutex tests |

## Common Errors

| Error or symptom | Cause | Correction |
|---|---|---|
| Multi-class events disappear | One class retained per cell | Keep unique classes separated by `;` |
| Wrong annotation over samples | Clinical rows used positionally | Match exact IDs and assert order |
| All-grey annotation | Clinical and mutation IDs do not match | Stop on missing matches; normalize IDs upstream deliberately |
| Misleading gene order | Default counts alteration events | Pass `row_order = order(-rowSums(mat != ''))` |
| Samples absent from maftools plot | They have no MAF row | Use a cohort-complete ComplexHeatmap matrix |
| All-empty matrix fails obscurely | No mapped event remains | Check classes/genes; the helper stops with context |
| Partial maftools annotation is grey | Incomplete `annotationColor` | Map every discrete level and numeric feature |
| comut top gene appears at bottom | Natural category order supplied | Reverse `category_order`; the script does this |
| comut continuous track errors | pandas 3 with comut 0.0.3 | Use pandas 2.x |
| TMB colors saturate | Hard-coded range | Use `(0, observed maximum)` |
| comut sample IDs overlap | Dense cohort exceeds 50 samples | Let the script hide them, or set `--sample-label-threshold` explicitly |
| comut accepts invalid burden values | TMB was not range-checked | Require finite, non-negative TMB; the script rejects invalid values |

## References

- Canisius S, Martens JWM, Wessels LFA. 2016. *Genome Biol* 17:261.
- Cerami E, Gao J, Dogrusoz U, et al. 2012. *Cancer Discov* 2:401-404.
- Gao J, Aksoy BA, Dogrusoz U, et al. 2013. *Sci Signal* 6:pl1.
- Gu Z, Eils R, Schlesner M. 2016. *Bioinformatics* 32:2847-2849.
- Mayakonda A, Lin DC, Assenov Y, Plass C, Koeffler HP. 2018. *Genome Res* 28:1747-1756.

## Related Skills

- data-visualization/heatmaps-clustering
- data-visualization/lollipop-protein-maps
- data-visualization/color-palettes
- clinical-databases/variant-prioritization
- variant-calling/variant-annotation
- copy-number/cnv-annotation
