---
name: bio-data-visualization-upset-plots
description: Build UpSet plots to visualize set intersections beyond 4 sets (where Venn fails) using ComplexUpset (modern, ggplot2-grammar) or the unmaintained UpSetR, with explicit cardinality vs degree sorting, attribute panels, and carefully scoped query highlighting. Use when comparing overlap across many gene sets, peak sets, variant lists, or any set membership matrix where Venn diagrams become illegible.
tool_type: mixed
primary_tool: ComplexUpset
license: MIT
author: GPTomics
---

# UpSet Plots

Use an UpSet plot for four or more sets. Each vertical bar is an **exclusive** membership combination; the left bars are set totals. Cardinality order answers “which combinations have most elements?” Degree order answers “which combinations involve the most sets?”

Tested stack: R 4.4.3, ComplexUpset 1.3.3, UpSetR 1.4.0, ggplot2 4.0.3; Python 3.12, `upsetplot==0.9.0`, `pandas>=2.2,<3`, numpy 2.5.3, matplotlib 3.11. `upsetplot 0.9.0` currently fails with pandas 3; create a pandas-2 environment rather than treating that error as a data problem. ComplexUpset 1.3.3 draws on ggplot2 4.0.3 in this stack; do not pin ggplot2 merely for this Skill.

## Input preflight (run before every R route)

This normalizes identifiers but deliberately preserves every named set. It warns when it removes missing/blank IDs or duplicates, and stops rather than silently dropping an empty set.

```r
prepare_sets <- function(sets) {
  stopifnot(is.list(sets), !is.null(names(sets)), all(nzchar(names(sets))))
  cleaned <- lapply(names(sets), function(nm) {
    x <- trimws(as.character(sets[[nm]]))
    bad <- is.na(x) | !nzchar(x)
    if (any(bad)) warning(sprintf("%s: removed %d NA/blank identifier(s)", nm, sum(bad)))
    x <- x[!bad]
    before <- length(x); x <- unique(x)
    if (length(x) != before) warning(sprintf("%s: removed %d duplicate identifier(s)", nm, before - length(x)))
    x
  })
  names(cleaned) <- names(sets)
  empty <- names(cleaned)[lengths(cleaned) == 0L]
  if (length(empty)) stop("empty set(s) after preflight: ", paste(empty, collapse = ", "),
                          ". Fix the input or omit them deliberately before plotting.")
  message("set sizes: ", paste(sprintf("%s=%d", names(cleaned), lengths(cleaned)), collapse = "; "),
          " | union=", length(unique(unlist(cleaned, use.names = FALSE))))
  cleaned
}

sets <- prepare_sets(sets)
```

For a pre-existing binary membership table, reject or resolve missing membership/ID values before plotting; do not let missing IDs create an apparent overlap. `fromList()` and the ComplexUpset `%in%` construction below both operate on de-duplicated sets, so duplicated IDs do **not** inflate counts. Python `from_contents()` rejects duplicate IDs instead.

## ComplexUpset (modern default)

The planted example has unequal exclusive intersections: `A` only=6, `B` only=4, `C` only=2, `A-B`=5, `A-C`=3, `B-C`=2, `A-B-C`=4. It makes sorting and annotation effects visible.

```r
library(ComplexUpset)
library(ggplot2)

sets <- list(
  SetA = c(paste0("g", 1:6), paste0("ab", 1:5), paste0("ac", 1:3), paste0("abc", 1:4)),
  SetB = c(paste0("b", 1:4), paste0("ab", 1:5), paste0("bc", 1:2), paste0("abc", 1:4)),
  SetC = c(paste0("c", 1:2), paste0("ac", 1:3), paste0("bc", 1:2), paste0("abc", 1:4))
)
sets <- prepare_sets(sets)
elements <- sort(unique(unlist(sets, use.names = FALSE)))
df <- data.frame(element = elements)
for (s in names(sets)) df[[s]] <- df$element %in% sets[[s]]

p <- ComplexUpset::upset(
  df, intersect = names(sets), n_intersections = 20,
  sort_intersections = "descending", sort_intersections_by = "cardinality",
  base_annotations = list("Intersection size" = intersection_size(counts = TRUE, text = list(size = 3))),
  themes = upset_modify_themes(list("Intersection size" = theme(panel.grid = element_blank())))
)
ggsave("upset_complex.pdf", p, width = 9, height = 5, device = cairo_pdf)
```

Use `sort_intersections_by="degree", sort_intersections="ascending"` for 1-set, then 2-set, then 3-set groups. To remove single-set combinations, use `min_degree=2` (or an explicit `intersections=list(...)`); `mode="intersect"` changes the count definition to inclusive set membership and is not a one-set filter. `n_intersections` limits the rendered rank, not the theoretical `2^N - 1` possibilities: 10 real sets can have far fewer non-empty combinations.

## Queries and attribute panels

Use a **single** `upset_query()` per figure with this ComplexUpset release. A single existing target renders correctly; multiple non-adjacent query geometries can span unrelated bars, so make separate figures for separate highlighted intersections. A query must target a non-empty exclusive combination in the supplied data.

```r
# SetA-SetB is present in the planted data; one highlight per figure.
p_query <- ComplexUpset::upset(
  df, intersect = names(sets),
  queries = list(upset_query(intersect = c("SetA", "SetB"), color = "#D55E00", fill = "#D55E00",
                             only_components = c("intersections_matrix", "Intersection size")))
)

set.seed(20260923)
df$log2FC <- rnorm(nrow(df), sd = 1.5)
df$significant <- df$log2FC > 1
p_attrs <- ComplexUpset::upset(
  df, intersect = names(sets),
  annotations = list(
    "log2FC" = ggplot(mapping = aes(x = intersection, y = log2FC)) + geom_boxplot() + theme_classic(),
    "Significant fraction" = ggplot(mapping = aes(x = intersection, fill = significant)) +
      geom_bar(position = "fill") +
      scale_fill_manual(values = c("TRUE" = "#D55E00", "FALSE" = "grey80")) + theme_classic()
  )
)
```

## UpSetR (legacy reproducibility only)

Always namespace it: both packages export `upset()`. Preserve all cleaned sets with `nsets=length(sets)`; a smaller value silently chooses only the largest sets.

```r
upset_data <- UpSetR::fromList(sets)
UpSetR::upset(upset_data, nsets = length(sets), nintersects = 20,
              order.by = "freq", decreasing = TRUE,
              mainbar.y.label = "Elements in exclusive intersection",
              sets.x.label = "Elements per set")
```

## Python: upsetplot 0.9.0

`show_counts=True` raises during export on the tested numpy 2.5 stack. Use `show_counts=False` and label the existing intersection bars yourself. `max_subset_rank`, not `intersection_plot_elements`, limits intersections by rank, but ties at the cutoff can retain more bars than the numeric rank. Pre-filter the intersection series when an exact bar count is required. `present` alone styles every superset containing those sets; combine it with `absent` to select one exact membership.

```python
import matplotlib as mpl
mpl.rcParams["pdf.fonttype"] = 42
import matplotlib.pyplot as plt
from upsetplot import from_contents, UpSet

sets = {"SetA": ["g1", "g2", "g3"], "SetB": ["g2", "g3", "g4"], "SetC": ["g3", "g5"]}
data = from_contents(sets)  # pandas DataFrame; duplicate identifiers raise ValueError
upset = UpSet(data, subset_size="count", show_counts=False, sort_by="cardinality",
              sort_categories_by="cardinality", max_subset_rank=20, facecolor="#0072B2")
# Exact SetA-and-SetB-only styling: state every other set is absent.
upset.style_subsets(present=["SetA", "SetB"], absent=["SetC"], facecolor="#D55E00")
fig = plt.figure(figsize=(8, 5))
axes = upset.plot(fig=fig)
for bar in axes["intersections"].patches:
    height = bar.get_height()
    if height:
        axes["intersections"].annotate(f"{height:g}", (bar.get_x() + bar.get_width() / 2, height),
                                        ha="center", va="bottom", fontsize=8)
fig.savefig("upset.pdf", bbox_inches="tight")
plt.close(fig)
```

## Checks before interpretation

- Print each cleaned set size and union size. In every plotted route, compare the sum of exclusive bar heights with the union size; it must agree when all non-empty intersections are rendered.
- Never claim that all `2^N-1` combinations are drawn: only non-empty combinations exist in the data, and rank/filter arguments further reduce them.
- Use `sort_by="cardinality"` / `sort_intersections_by="cardinality"` for largest first. Degree order in upsetplot is ascending by default; use `sort_by="-degree"` only when high-degree combinations must come first.
- Keep 15–25 rendered combinations for a readable static figure; use `min_degree`, `min_subset_size`, or an explicit intersection list to address the scientific question.

## Common errors

| Symptom | Cause | Correction |
|---|---|---|
| Missing/blank IDs appear shared | Identifier hygiene was skipped | Run `prepare_sets()` and repair the source IDs. |
| An expected set vanished | It was empty, or UpSetR received a too-small `nsets` | Stop on empty input; use `nsets=length(sets)`. |
| Two ComplexUpset highlights span other bars | Multi-query geometry limitation in 1.3.3 | Render one existing query per figure. |
| Python `Invalid RGBA argument: nan` | pandas 3 with upsetplot 0.9.0 | Use pandas >=2.2,<3 or upgrade only after testing a newer upsetplot. |
| Python highlight colors supersets too | `present` does not mean exact membership | Add `absent` for every other set. |

## References

- Lex A, Gehlenborg N, Strobelt H, Vuillemot R, Pfister H. 2014. UpSet. *IEEE TVCG* 20(12):1983-1993.
- Conway JR, Lex A, Gehlenborg N. 2017. UpSetR. *Bioinformatics* 33(18):2938-2940.
- Krassowski M. ComplexUpset. https://github.com/krassowski/complex-upset

## Related Skills

- data-visualization/heatmaps-clustering
- pathway-analysis/go-enrichment
- differential-expression/de-results
