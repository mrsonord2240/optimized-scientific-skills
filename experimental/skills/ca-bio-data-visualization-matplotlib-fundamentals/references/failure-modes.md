# matplotlib failure modes

## Default Type-3 fonts rejected by journals

**Trigger:** Default `pdf.fonttype=3` (PostScript Type 3 glyphs as drawing operators).

**Mechanism:** Type-3 glyphs are not searchable or selectable; many journals reject.

**Symptom:** Submission rejected at automated check; "Type 3 fonts not permitted."

**Fix:** `mpl.rcParams['pdf.fonttype']=42` AND `ps.fonttype=42`. Verify with `pdffonts figure.pdf` showing `TrueType`.

## tight_layout does not follow later changes

**Trigger:** `fig.tight_layout()` called early, then a label is lengthened or the figure is resized; or `fig.colorbar(im, ax=[a, b])` spanning several axes.

**Mechanism:** `tight_layout()` is a one-shot adjustment of the subplot margins at call time. A layout engine re-runs on every draw.

**Symptom:** Measured on matplotlib 3.11.2 (89 x 60 mm): axis labels clipped by the page edge after the later change; multi-axes colorbars warn "Axes that are not compatible with tight_layout". A single colorbar on one axes did not clip.

**Fix:** `plt.subplots(layout='constrained')` at creation. `fig.set_layout_engine('constrained')` also works, but only before any colorbar exists; afterwards it raises `ZeroDivisionError` on 3.11.2.

## pyplot state-machine in multi-subplot

**Trigger:** `plt.xlabel(...)` after `plt.subplots(2, 3)`.

**Mechanism:** pyplot calls modify the *current* axes — usually the last created. Multi-subplot code becomes order-dependent.

**Symptom:** Wrong subplot gets the label.

**Fix:** Use `ax.set_xlabel(...)` with explicit axes reference.

## Vector scatter at large N bloats the PDF

**Trigger:** Vector scatter at large N.

**Mechanism:** Each scatter point is a vector path.

**Symptom:** Measured on matplotlib 3.11.2 at 89 mm: 100000 points = 1.5 MB vector vs 24 KB rasterized (63x); slow to open in Illustrator and viewers.

**Fix:** `rasterized=True` on the scatter call. Keep axes and text vector.

## seaborn FacetGrid vs Axes return-type confusion

**Trigger:** `g = sns.displot(...)`; calling `g.set_xlabel('x')` fails.

**Mechanism:** displot returns FacetGrid; needs `.set_axis_labels(x, y)` or per-axes iteration.

**Symptom:** AttributeError on .set_xlabel.

**Fix:** Use `set_axis_labels` for FacetGrid; `set_xlabel` for Axes. Switch to axes-level `sns.histplot(ax=ax)` to get Axes-API behavior.

## figsize in inches when mm was intended

**Trigger:** `figsize=(89, 70)` thinking mm; matplotlib expects inches.

**Mechanism:** Default figure unit is inches.

**Symptom:** Figure is 89 inches wide.

**Fix:** Convert: `figsize=(89/25.4, 70/25.4)` for mm input.

## Colorbar over-fills the axes

**Trigger:** Default `fig.colorbar(im, ax=ax)`.

**Mechanism:** Colorbar takes the same height as the axes; on small subplots dominates.

**Symptom:** Subplot looks squished.

**Fix:** `fig.colorbar(im, ax=ax, shrink=0.6, aspect=20)`; or use `make_axes_locatable` for fine control.

## Vector grid + rasterized scatter mixed properly

**Trigger:** Want vector axes + raster scatter; save as PDF.

**Mechanism:** Rasterization is per artist; without `rasterized=True` everything stays vector.

**Symptom:** Whole plot rasterized; axis text blurry on zoom.

**Fix:** Per-element `rasterized=True` on scatter only; axes and text stay vector. To rasterize every artist below a zorder, call `ax.set_rasterization_zorder(1)` (an Axes method; artists with `zorder < 1` are rasterized).

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| PDF rejected by journal | Type-3 fonts | `pdf.fonttype=42` |
| Subplots overlap | No constrained layout | `plt.subplots(layout='constrained')` |
| Wrong subplot labeled | pyplot state-machine | Use ax.set_xlabel explicitly |
| Multi-MB PDF | Vector scatter at large N | `rasterized=True` on scatter |
| Figure too big | mm interpreted as inches | Divide by 25.4 |
| Colorbar dominates | Default size | `shrink=0.6, aspect=20` |
| seaborn .set_xlabel fails | FacetGrid not Axes | `g.set_axis_labels(x, y)` |
| Journal-width page is ~92 mm | `bbox_inches='tight'` re-crops the page | Drop it; rely on `layout='constrained'` |
| SVG text not editable | `svg.fonttype='path'` (default) | `svg.fonttype='none'` |
| Axes spine missing | Wrong API | `ax.spines[['top','right']].set_visible(False)` |
