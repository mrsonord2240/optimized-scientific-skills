# ggplot2 failure modes

## Default ggsave fonts not embedded

**Trigger:** `ggsave('out.pdf', p)` without `device = cairo_pdf`.

**Mechanism:** The default pdf() device does not embed Helvetica/Symbol (verified with `pdffonts`); cairo_pdf embeds TrueType.

**Symptom:** Reviewer or coauthor opens PDF; text renders in wrong font; journal rejects.

**Fix:** Always `device = cairo_pdf` for PDF saves.

## Mapping vs constant aesthetic confusion

**Trigger:** `geom_point(aes(color = 'red'))` — string 'red' becomes a categorical mapping.

**Mechanism:** `aes()` interprets its arguments as variables; 'red' becomes a 1-level factor and gets mapped to the first default hue (#F8766D, salmon), not red.

**Symptom:** Points appear salmon (not red) with a legend showing "red" as a category.

**Fix:** Move outside aes: `geom_point(color = 'red')` for a constant; keep inside for a mapping.

## linewidth vs size for lines

**Trigger:** `geom_line(size = 0.5)` in ggplot2 3.4+.

**Mechanism:** ggplot2 3.4+ renamed line-width control from `size` to `linewidth`; `size` still works for points.

**Symptom:** Warning "Using `size` aesthetic for lines was deprecated"; lines render but warning.

**Fix:** `geom_line(linewidth = 0.5)`. `geom_point(size = 1)` is correct.

## facet_wrap scales = 'free' confuses cross-panel comparison

**Trigger:** `facet_wrap(~ var, scales = 'free')` for figures intended to compare across panels.

**Mechanism:** Each panel has its own scale; visual comparison invalid.

**Symptom:** Reviewer asks "why are these heights different?"

**Fix:** Use `scales = 'fixed'` (default) when cross-panel comparison matters; use `'free_y'` only when panels are inherently different scales.

## aes_string deprecated

**Trigger:** `aes_string(x = 'PC1', y = 'PC2')` for programmatic plotting.

**Mechanism:** Deprecated since ggplot2 3.0; emits warning.

**Symptom:** Deprecation warning in script log.

**Fix:** `aes(x = .data[['PC1']], y = .data[['PC2']])` OR `aes(x = !!sym(x_var))`.

## ggrepel max.overlaps default drops labels

**Trigger:** `geom_text_repel(aes(label = label))` where a label overlaps more than 10 other labels or points (dense volcano labels, long IDs).

**Mechanism:** Default `max.overlaps = 10`; a label overlapping more items than that is dropped with no warning or message (ggrepel 0.9.8 reports it only when `verbose = TRUE`). Few labels in a sparse plot are unaffected.

**Symptom:** Some labeled genes are missing and nothing in the log says so (measured: 120 of 300 labels drawn at the default, 300 with `max.overlaps = Inf`).

**Fix:** `geom_text_repel(aes(label = label), max.overlaps = Inf, seed = 1)` OR `options(ggrepel.max.overlaps = Inf)` at script top; then count drawn labels and open the figure for overlaps.

## Saving with size in inches but intended mm

**Trigger:** `ggsave('out.pdf', p, width = 89, height = 70)` thinking mm.

**Mechanism:** Default `units = 'in'`.

**Symptom:** Figure is 89 inches wide — too large to open in Illustrator.

**Fix:** `units = 'mm'` explicit. Nature single column = 89mm; double column = 183mm.
