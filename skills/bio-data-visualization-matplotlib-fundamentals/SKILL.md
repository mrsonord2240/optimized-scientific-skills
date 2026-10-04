---
name: bio-data-visualization-matplotlib-fundamentals
description: Use when producing publication figures in Python with matplotlib, such as RNA-seq scatter plots, single-cell embeddings, or general biological plots.
tool_type: python
primary_tool: matplotlib
goal_approach_exempt: true
license: MIT
author: GPTomics
category: Data Analysis
---

## Version Compatibility

Reference examples tested with: matplotlib 3.8+ (re-checked on matplotlib 3.11.2), seaborn 0.13+ (0.13.2), numpy 1.26+ (2.5.3), pandas 2.2+ (3.0.6), cmcrameri 1.9+ (1.10).

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed package and adapt the example to match the actual API rather than retrying. Renamed keywords seen on 3.11: `Axes.boxplot(labels=)` is now `tick_labels=`.

# matplotlib Fundamentals

**"Make a publication figure in Python"** -> Build via the **object-oriented Figure/Axes API** (not pyplot state-machine), with constrained layout for axes alignment, `pdf.fonttype=42` for journal-compliant TrueType fonts, CVD-safe palettes, and rasterized point layers for large scatter. The pyplot interface is for notebook scratch; the Figure/Axes API is for reproducible figures.

- Python: `fig, ax = plt.subplots()` -> `ax.scatter` / `ax.plot` / `ax.bar`; `seaborn.objects` (new grammar API) for ggplot-like

## The Three Modern Defaults

1. **Object-oriented API** — `fig, ax = plt.subplots(figsize=(4, 3))` then `ax.scatter(x, y)`, `ax.set_xlabel(...)`. The pyplot state-machine (`plt.scatter`, `plt.xlabel`) hides which axes are being modified and breaks in multi-subplot figures.

2. **Constrained layout** — `plt.subplots(layout='constrained')` automatically prevents axis-label clipping and tight-packs subplots. It replaces the older `tight_layout()`, is opt-in (the default layout engine is none), and supersedes the `constrained_layout=True` keyword.

3. **Type-42 (TrueType) font embedding** — `plt.rcParams['pdf.fonttype']=42` produces searchable/editable PDF text. Default Type-3 PostScript glyphs are not searchable and **rejected by Nature, IEEE, ACM, and many other publishers**.

## Standard Setup for Publication

```python
import matplotlib.pyplot as plt
import matplotlib as mpl

# rcParams for publication compliance
mpl.rcParams.update({
    'pdf.fonttype': 42,                 # TrueType -- searchable PDFs
    'ps.fonttype': 42,                  # TrueType in EPS
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 7,                     # Nature requires 5-7 pt body text
    'axes.labelsize': 7,
    'axes.titlesize': 7,
    'xtick.labelsize': 6,
    'ytick.labelsize': 6,
    'legend.fontsize': 6,
    'figure.dpi': 100,                  # display
    'savefig.dpi': 300,                 # save
    'axes.linewidth': 0.5,
    'xtick.major.width': 0.5,
    'ytick.major.width': 0.5,
    'lines.linewidth': 1.0,
    'patch.linewidth': 0.5,
})
```

## Figure / Axes API

```python
import matplotlib.pyplot as plt

# Single axes
fig, ax = plt.subplots(figsize=(89/25.4, 70/25.4),       # 89mm x 70mm in inches; Nature single col
                       layout='constrained')
ax.scatter(x, y, c='#0072B2', s=10, alpha=0.7, edgecolors='none', rasterized=True)
ax.set_xlabel('PC1 (45%)')
ax.set_ylabel('PC2 (12%)')
ax.spines[['top', 'right']].set_visible(False)
fig.savefig('scatter.pdf')

# Grid of axes
fig, axes = plt.subplots(2, 3, figsize=(180/25.4, 100/25.4),  # 180mm double col
                          layout='constrained')
for ax, (label, panel_data) in zip(axes.flat, panel_data.items()):
    ax.plot(panel_data['x'], panel_data['y'])
    ax.text(-0.15, 1.05, label, transform=ax.transAxes, fontsize=7, fontweight='bold', va='top')  # panel tag
```

## seaborn Integration

```python
import seaborn as sns

# seaborn shares the matplotlib Figure/Axes -- pass ax= argument
fig, ax = plt.subplots(figsize=(4, 3), layout='constrained')
# Bind colours to category names: a list palette follows data order and can swap Up and Down.
sig_palette = {'NS': '#999999', 'Down': '#0072B2', 'Up': '#D55E00'}
sns.scatterplot(data=df, x='log_fold_change', y='neg_log_p',
                hue='significance', hue_order=['NS', 'Down', 'Up'], palette=sig_palette,
                s=10, alpha=0.7, ax=ax, rasterized=True)

# seaborn 0.13+ has the `objects` grammar interface (ggplot-like)
import seaborn.objects as so
(so.Plot(df, x='log_fold_change', y='neg_log_p')
   .add(so.Dots(pointsize=2), color='significance')
   .scale(color=sig_palette))     # dict, not list
# so.Plot restyles with its own theme; pass .theme({**sns.axes_style('ticks'), **mpl.rcParams}) to keep the rcParams above.
```

**Return-type gotcha:** seaborn axes-level functions (`scatterplot`, `boxplot`, `barplot`) return Axes. Figure-level (`displot`, `relplot`, `catplot`) return FacetGrid — needs `.set_axis_labels(x, y)` not `.set_xlabel(x)`.

## Color and Palette

```python
# CVD-safe categorical
okabe_ito = ['#E69F00', '#56B4E9', '#009E73', '#F0E442',
             '#0072B2', '#D55E00', '#CC79A7', '#000000']

# Perceptually-uniform sequential (Crameri batlow / viridis cividis)
import numpy as np
from cmcrameri import cm as cmc
fig, ax = plt.subplots(layout='constrained')
im = ax.imshow(data, cmap=cmc.batlow)                    # or cmap='viridis' (built-in)
fig.colorbar(im, ax=ax, shrink=0.6, aspect=20)

# Diverging symmetric for LFC / z-score
vmax = np.quantile(np.abs(data), 0.99)
im = ax.imshow(data, cmap='RdBu_r', vmin=-vmax, vmax=vmax)   # symmetric
```

See `data-visualization/color-palettes` for full palette decision tree.

## Saving

```python
# PDF for vector text + raster scatter (best of both)
fig.savefig('figure.pdf')

# PNG for raster (web, presentations); dpi comes from rcParams savefig.dpi
fig.savefig('figure.png', dpi=300)

# TIFF for some journals
fig.savefig('figure.tiff', dpi=300, pil_kwargs={'compression': 'tiff_lzw'})

# SVG: text is converted to paths by default; keep it editable as <text> (the viewer needs the font installed)
with mpl.rc_context({'svg.fonttype': 'none'}):
    fig.savefig('figure.svg')
```

## Guardrails

- Set `pdf.fonttype=42` and `ps.fonttype=42` before any PDF/EPS save; verify with `pdffonts figure.pdf` (no Type 3 rows).
- Use `layout='constrained'`, not `tight_layout()`, when colorbars or shared axes are present.
- Do not set `bbox_inches='tight'` (or `savefig.bbox`) for journal-width figures: it re-crops the page, so 89 mm becomes ~92 mm. Check the saved page size in mm.
- Map categories to colours with a dict palette (`palette={'Up': ...}`, plus `hue_order`), never a list in assumed data order.
- Address axes explicitly (`ax.set_xlabel`), never `plt.xlabel` after `plt.subplots(n, m)`.
- `rasterized=True` on scatter/imshow above ~1000 elements; axes and text stay vector.
- `figsize` is inches: `mm / 25.4`. Colorbar: `shrink=0.6, aspect=20` on small axes.
- seaborn figure-level functions (`displot`, `relplot`, `catplot`) return a FacetGrid (`.set_axis_labels`); axes-level ones accept `ax=`.

Trigger/mechanism/symptom/fix and the error table: `references/failure-modes.md`. Chart-type and axis-formatting recipes: `references/chart-recipes.md`.

## Resources

- `scripts/matplotlib_phd.py` — runnable end-to-end example (rcParams, 89 mm scatter, 2x3 grid, Crameri heatmap, seaborn, seaborn.objects; its legend recipe uses the private `Plot.plot()._figure` (seaborn 0.13.2 only) and a fixed 22% page reserve that a longer legend title overflows, a 29-character title covered 13 mm of data; the public alternative is `Plot.on(fig)`, move `fig.legends[0]`, `fig.subplots_adjust(right=0.76)`). Writes five PDFs at their stated mm sizes to the working directory; needs cmcrameri.
- `usage-guide.md` — prerequisites and example prompts.

## References

- Hunter JD. 2007. Matplotlib: A 2D graphics environment. *Comput Sci Eng* 9(3):90-95.
- Rougier NP, Droettboom M, Bourne PE. 2014. Ten simple rules for better figures. *PLOS Comp Biol* 10(9):e1003833.
- Waskom ML. 2021. seaborn: statistical data visualization. *J Open Source Softw* 6(60):3021.

## Related Skills

- data-visualization/color-palettes - Palette selection
- data-visualization/multipanel-figures - GridSpec and patchwork-equivalent layouts
- data-visualization/distribution-plots - seaborn boxplot/violin/raincloud
- data-visualization/heatmaps-clustering - seaborn.clustermap
- data-visualization/volcano-and-ma-plots - matplotlib scatter for volcano
- reporting/figure-export - DPI / format / journal-spec details
